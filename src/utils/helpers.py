#!/usr/bin/env python3
"""Helper utilities - OPTIMIZED VERSION: Faster processing, no partial successes"""
import gzip
import shutil
import hashlib
import subprocess
import boto3
import logging
import re
import mimetypes
import psutil
import os
from pathlib import Path
from Bio import SeqIO
from botocore.config import Config
from boto3.s3.transfer import TransferConfig
from botocore.exceptions import ClientError, EndpointConnectionError
from concurrent.futures import ProcessPoolExecutor, as_completed

logger = logging.getLogger(__name__)


# =========================
# File validation & parsing
# =========================
class FileValidator:
    @staticmethod
    def validate_and_fix(file_path: Path):
        if not file_path.exists() or file_path.stat().st_size < 100:
            return None
        # if gz, decompress for validation
        if file_path.suffix == ".gz":
            decompressed = file_path.with_suffix("")
            if not decompressed.exists():
                try:
                    with gzip.open(file_path, "rb") as f_in, open(decompressed, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                except Exception:
                    return None
            file_path = decompressed
        try:
            cnt = sum(1 for _ in SeqIO.parse(file_path, "fasta"))
            if cnt > 0:
                return file_path
        except Exception:
            pass
        return None

    @staticmethod
    def _validate_gzip_ok(gz_path: Path) -> bool:
        try:
            with gzip.open(gz_path, "rb") as fh:
                _ = fh.read(1024)
            return True
        except Exception as e:
            logger.error(f"Gzip validation failed for {gz_path}: {e}")
            return False

    @staticmethod
    def ensure_gz(path: Path) -> Path | None:
        """Ensure file is gzipped, return gzipped path or None if failed"""
        if not path or not path.exists():
            return None
            
        if path.suffix == ".gz":
            return path if FileValidator._validate_gzip_ok(path) else None
        
        gz_path = path.with_suffix(path.suffix + ".gz")
        
        # Check if already compressed and valid
        if (
            gz_path.exists()
            and gz_path.stat().st_size > 0
            and gz_path.stat().st_mtime >= path.stat().st_mtime
            and FileValidator._validate_gzip_ok(gz_path)
        ):
            return gz_path
        
        # Compress the file
        try:
            with open(path, "rb") as fin, gzip.open(gz_path, "wb") as fout:
                shutil.copyfileobj(fin, fout)
        except Exception as e:
            logger.error(f"Compression failed for {path} -> {gz_path}: {e}")
            return None
        
        # Validate compressed file
        if FileValidator._validate_gzip_ok(gz_path):
            # Delete original uncompressed file after successful compression
            try:
                path.unlink()
                logger.info(f"Deleted original uncompressed file: {path.name}")
            except Exception as e:
                logger.warning(f"Failed to delete original file {path}: {e}")
            return gz_path
        else:
            return None

    @staticmethod
    def parse_species_family(file_path: Path):
        try:
            for rec in SeqIO.parse(file_path, "fasta"):
                desc = rec.description
                if "[" in desc and "]" in desc:
                    return None, None
                parts = desc.split(maxsplit=1)
                if len(parts) < 2:
                    break
                tail = parts[1]
                tail = re.split(
                    r",|\bchromosome\b|\bscaffold\b|\bcontig\b|\bstrain\b|\bcomplete\b",
                    tail, flags=re.I
                )[0].strip()
                toks = tail.split()
                if len(toks) >= 2:
                    species = f"{toks[0]} {toks[1]}"
                elif toks:
                    species = toks[0]
                else:
                    species = None
                if not species or "[" in species or "]" in species or len(species.split()) < 2:
                    return None, None
                return species, None
        except Exception:
            pass
        return None, None

    @staticmethod
    def compute_sequence_hash(file_path: Path):
        try:
            hashes = []
            for rec in SeqIO.parse(file_path, "fasta"):
                seq = str(rec.seq).upper().replace("N", "")
                if seq:
                    hashes.append(hashlib.sha256(seq.encode()).hexdigest())
            if hashes:
                hashes.sort()
                return hashlib.sha256("|".join(hashes).encode()).hexdigest()
        except Exception:
            pass
        return None


# =========================
# Disk management
# =========================
class DiskManager:
    @staticmethod
    def bytes_free(path: str = "/") -> int:
        try:
            return psutil.disk_usage(path).free
        except Exception:
            st = os.statvfs(path)
            return st.f_bavail * st.f_frsize

    @staticmethod
    def ensure_free_space(min_free_gb: int, purge_paths: list[Path]) -> bool:
        free_gb = DiskManager.bytes_free("/") / (1024**3)
        if free_gb >= min_free_gb:
            return True
        # purge splits and .part under purge_paths
        for root in purge_paths:
            try:
                if root.is_dir():
                    for sp in root.rglob("splits"):
                        shutil.rmtree(sp, ignore_errors=True)
                    for pr in root.rglob("*.part"):
                        try:
                            pr.unlink()
                        except:
                            pass
            except:
                pass
        free_gb2 = DiskManager.bytes_free("/") / (1024**3)
        return free_gb2 >= min_free_gb

    @staticmethod
    def cleanup_files(files):
        """Clean up a list of files and their parent directories if empty"""
        for f in files:
            try:
                if f.exists():
                    f.unlink()
                parent = f.parent
                if parent.is_dir() and not any(parent.iterdir()):
                    parent.rmdir()
            except Exception as e:
                logger.debug(f"Failed to cleanup {f}: {e}")

    @staticmethod
    def cleanup_directory(directory: Path):
        """Remove directory and all contents"""
        if directory.exists():
            try:
                shutil.rmtree(directory)
                logger.info(f"Cleaned up directory: {directory}")
            except Exception as e:
                logger.warning(f"Failed to cleanup directory {directory}: {e}")


# =========================
# S3 Uploads (hardened)
# =========================
class S3Manager:
    def __init__(self, config):
        s3_cfg = config["aws"]["s3"]
        self.bucket = s3_cfg["bucket_name"]
        self.final_prefix = s3_cfg["final_prefix"].rstrip("/") + "/"
        self.proteins_prefix = s3_cfg["proteins_prefix"].rstrip("/") + "/"

        region = config["aws"].get("region", "us-east-1")

        self.s3 = boto3.client(
            "s3",
            region_name=region,
            config=Config(
                s3={"addressing_style": "virtual"},
                retries={"max_attempts": 8, "mode": "standard"},
                signature_version="s3v4",
                connect_timeout=10,
                read_timeout=300,
            ),
        )

        self.tcfg = TransferConfig(
            multipart_threshold=500 * 1024 * 1024,
            multipart_chunksize=64 * 1024 * 1024,
            max_concurrency=4,
            use_threads=True,
        )

    def _guess_meta(self, path: Path):
        ctype, _ = mimetypes.guess_type(str(path))
        extra = {}
        if path.suffix == ".gz":
            # DO NOT set ContentEncoding for .gz payloads — serve as binary gzip
            ctype = "application/gzip"
        extra["ContentType"] = ctype or "application/octet-stream"
        return extra

    def _upload(self, local_path: Path, s3_key: str) -> str | None:
        """Upload file to S3 and return S3 URI on success"""
        if not local_path or not local_path.exists() or local_path.stat().st_size == 0:
            logger.error(f"S3 upload skipped: local file missing/empty: {local_path}")
            return None
        
        extra_args = self._guess_meta(local_path)
        try:
            self.s3.upload_file(
                Filename=str(local_path),
                Bucket=self.bucket,
                Key=s3_key,
                ExtraArgs=extra_args,
                Config=self.tcfg,
            )
            
            # Verify upload
            head = self.s3.head_object(Bucket=self.bucket, Key=s3_key)
            etag = head.get("ETag", "").strip('"')
            size = head.get("ContentLength", 0)
            
            # Cheap gzip magic check
            if str(local_path).endswith(".gz"):
                try:
                    obj = self.s3.get_object(Bucket=self.bucket, Key=s3_key, Range="bytes=0-1")
                    magic = obj["Body"].read(2)
                    if magic != b"\x1f\x8b":
                        logger.warning(f"S3 object {s3_key} not gzip magic; got {magic!r}")
                except Exception as e:
                    logger.warning(f"Gzip magic check skipped for {s3_key}: {e}")
            
            logger.info(f"S3 ✓ {self.bucket}/{s3_key}  (ETag={etag}, size={size:,} bytes)")
            return f"s3://{self.bucket}/{s3_key}"
            
        except (ClientError, EndpointConnectionError) as e:
            logger.error(f"S3 upload failed for {s3_key}: {e}", exc_info=True)
            return None
        except Exception as e:
            logger.error(f"S3 unexpected error for {s3_key}: {e}", exc_info=True)
            return None

    def upload_genome(self, local_path: Path, source: str, accession: str) -> str | None:
        """Upload genome file to S3"""
        s3_key = f"{self.final_prefix}{source}/genomes/{accession}/{local_path.name}"
        return self._upload(local_path, s3_key)

    def upload_proteins(self, local_path: Path, source: str, accession: str) -> str | None:
        """Upload protein file to S3"""
        s3_key = f"{self.proteins_prefix}{source}/{accession}/{local_path.name}"
        return self._upload(local_path, s3_key)


# =========================
# Prodigal & metagenome helpers - OPTIMIZED VERSION
# =========================
class ProdigalRunner:
    @staticmethod
    def predict_proteins(genome_path: Path, is_metagenome: bool = False):
        """Run Prodigal on a genome file"""
        protein_path = genome_path.with_suffix(".faa")
        cmd = ["prodigal", "-i", str(genome_path), "-a", str(protein_path), "-o", "/dev/null", "-q"]
        if is_metagenome:
            cmd.extend(["-p", "meta"])
        
        try:
            subprocess.run(cmd, check=True, capture_output=True, timeout=300)
            if not protein_path.exists():
                return None
            
            # Compress and delete original
            gz = FileValidator.ensure_gz(protein_path)
            if gz:
                return gz
            else:
                logger.error(f"Protein compression failed for {protein_path}")
                return None
                
        except Exception as e:
            logger.error(f"Prodigal failed for {genome_path}: {e}")
            return None

    @staticmethod
    def _run_prodigal_on_split(split_path: Path):
        """Worker function for parallel Prodigal execution"""
        try:
            return ProdigalRunner.predict_proteins(split_path, is_metagenome=True)
        except Exception as e:
            logger.error(f"Prodigal worker failed for {split_path.name}: {e}")
            return None

    @staticmethod
    def process_metagenome(fasta_path: Path, output_dir: Path):
        """
        OPTIMIZED metagenome processing:
        1. Pre-filter short sequences (< 200bp - unlikely to have complete ORFs)
        2. Use SIZE-based splitting for even distribution
        3. Parallel Prodigal execution with ProcessPoolExecutor
        4. FAIL FAST - no partial successes
        5. Faster due to better parallelization
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        split_dir = output_dir / "splits"
        split_dir.mkdir(exist_ok=True)
        
        try:
            # Step 1: Pre-filter and get statistics
            logger.info(f"Pre-processing {fasta_path.name}...")
            filtered_path = output_dir / f"{fasta_path.stem}_filtered.fasta"
            
            total_seqs = 0
            total_bp = 0
            kept_seqs = 0
            kept_bp = 0
            filtered_records = []
            
            MIN_CONTIG_LENGTH = 200  # Minimum length for potential ORFs
            
            for rec in SeqIO.parse(fasta_path, 'fasta'):
                total_seqs += 1
                total_bp += len(rec.seq)
                if len(rec.seq) >= MIN_CONTIG_LENGTH:
                    filtered_records.append(rec)
                    kept_seqs += 1
                    kept_bp += len(rec.seq)
            
            logger.info(f"Filtered: {kept_seqs}/{total_seqs} seqs ({kept_bp:,}/{total_bp:,} bp)")
            
            # Fail if insufficient data after filtering
            if kept_seqs < 10:
                logger.error(f"Too few sequences after filtering ({kept_seqs} seqs). Minimum: 10")
                shutil.rmtree(split_dir, ignore_errors=True)
                return None
            
            if kept_bp < 50000:  # Less than 50kb total
                logger.error(f"Insufficient bases after filtering ({kept_bp:,} bp). Minimum: 50kb")
                shutil.rmtree(split_dir, ignore_errors=True)
                return None
            
            # Write filtered sequences
            SeqIO.write(filtered_records, filtered_path, 'fasta')
            
            # Step 2: Determine splitting strategy
            # Small files: process directly without splitting
            if kept_seqs < 1000 or kept_bp < 500000:  # < 1000 seqs or < 500kb
                logger.info(f"Small metagenome, processing without splitting")
                protein_path = ProdigalRunner.predict_proteins(filtered_path, is_metagenome=True)
                shutil.rmtree(split_dir, ignore_errors=True)
                if filtered_path.exists():
                    filtered_path.unlink()
                return protein_path
            
            # Step 3: Size-based splitting (ensures even distribution)
            # Target: ~100kb per split, max 8 splits
            target_size_kb = 100
            num_splits = min(8, max(2, kept_bp // (target_size_kb * 1000)))
            chunk_size_kb = kept_bp // (num_splits * 1000)
            
            logger.info(f"Splitting into {num_splits} chunks (~{chunk_size_kb}kb each)")
            
            # Use seqkit split by size (more even than by parts)
            try:
                subprocess.run(
                    ["seqkit", "split2", "-s", f"{chunk_size_kb}k", "-O", str(split_dir), str(filtered_path)],
                    check=True,
                    capture_output=True,
                    timeout=300
                )
            except subprocess.CalledProcessError:
                # Fallback to split by parts if split2 not available
                logger.warning("seqkit split2 failed, falling back to split by parts")
                subprocess.run(
                    ["seqkit", "split", "-p", str(num_splits), "-O", str(split_dir), str(filtered_path)],
                    check=True,
                    timeout=300
                )
            
            split_files = sorted(split_dir.glob("*.fasta"))
            if not split_files:
                logger.error(f"No split files created for {fasta_path.name}")
                shutil.rmtree(split_dir, ignore_errors=True)
                if filtered_path.exists():
                    filtered_path.unlink()
                return None
            
            logger.info(f"Created {len(split_files)} split files")
            
            # Step 4: Validate all splits BEFORE running Prodigal
            for split_file in split_files:
                split_seqs = sum(1 for _ in SeqIO.parse(split_file, 'fasta'))
                split_bp = sum(len(r.seq) for r in SeqIO.parse(split_file, 'fasta'))
                
                if split_seqs < 10:
                    logger.error(f"Split {split_file.name} too small: {split_seqs} seqs")
                    shutil.rmtree(split_dir, ignore_errors=True)
                    if filtered_path.exists():
                        filtered_path.unlink()
                    return None
                
                if split_bp < 10000:  # Less than 10kb
                    logger.error(f"Split {split_file.name} insufficient data: {split_bp:,} bp")
                    shutil.rmtree(split_dir, ignore_errors=True)
                    if filtered_path.exists():
                        filtered_path.unlink()
                    return None
            
            # Step 5: Parallel Prodigal execution
            logger.info(f"Running Prodigal on {len(split_files)} splits in parallel...")
            
            protein_files = []
            with ProcessPoolExecutor(max_workers=min(8, len(split_files))) as executor:
                future_to_split = {
                    executor.submit(ProdigalRunner._run_prodigal_on_split, split_file): split_file 
                    for split_file in split_files
                }
                
                for future in as_completed(future_to_split):
                    split_file = future_to_split[future]
                    try:
                        result = future.result()
                        if result is None or not result.exists():
                            # FAIL FAST - any failure means total failure
                            logger.error(f"Prodigal failed on split {split_file.name}")
                            executor.shutdown(wait=False, cancel_futures=True)
                            shutil.rmtree(split_dir, ignore_errors=True)
                            if filtered_path.exists():
                                filtered_path.unlink()
                            return None
                        protein_files.append(result)
                    except Exception as e:
                        logger.error(f"Exception processing {split_file.name}: {e}")
                        executor.shutdown(wait=False, cancel_futures=True)
                        shutil.rmtree(split_dir, ignore_errors=True)
                        if filtered_path.exists():
                            filtered_path.unlink()
                        return None
            
            # Step 6: Merge all proteins (all must succeed to reach here)
            logger.info(f"Merging {len(protein_files)} protein files...")
            final_proteins = output_dir / f"{fasta_path.stem}_proteins.faa"
            
            with open(final_proteins, "w") as outf:
                for faa in sorted(protein_files):
                    if faa.suffix == '.gz':
                        with gzip.open(faa, 'rt') as inf:
                            outf.write(inf.read())
                    else:
                        with open(faa) as inf:
                            outf.write(inf.read())
            
            # Step 7: Cleanup and compress
            shutil.rmtree(split_dir, ignore_errors=True)
            if filtered_path.exists():
                filtered_path.unlink()
            
            if not final_proteins.exists() or final_proteins.stat().st_size == 0:
                logger.error(f"Merged protein file is empty")
                return None
            
            logger.info(f"✓ Successfully processed all {len(protein_files)} splits")
            
            # Compress and return
            gz = FileValidator.ensure_gz(final_proteins)
            if gz:
                return gz
            else:
                logger.error(f"Protein compression failed")
                return None
                
        except Exception as e:
            logger.error(f"Metagenome processing failed for {fasta_path}: {e}")
            if split_dir.exists():
                shutil.rmtree(split_dir, ignore_errors=True)
            return None
