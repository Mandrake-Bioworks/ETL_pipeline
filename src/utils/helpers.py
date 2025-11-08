#!/usr/bin/env python3
"""Helper utilities for file handling, validation, disk management, S3 uploads, and Prodigal"""
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
# Prodigal & metagenome helpers
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
    def process_metagenome(fasta_path: Path, output_dir: Path):
        """Process metagenome: split, run Prodigal on splits, merge results"""
        output_dir.mkdir(parents=True, exist_ok=True)
        split_dir = output_dir / "splits"
        split_dir.mkdir(exist_ok=True)
        
        try:
            # Split into chunks
            subprocess.run(
                ["seqkit", "split", "-p", "8", "-O", str(split_dir), str(fasta_path)],
                check=True, 
                timeout=600
            )
            
            split_files = list(split_dir.glob("*.fasta"))
            if not split_files:
                logger.warning(f"No split files created for {fasta_path}")
                return None
            
            # Run Prodigal on each split
            for split_file in split_files:
                ProdigalRunner.predict_proteins(split_file, is_metagenome=True)
            
            # Merge all protein predictions
            final_proteins = output_dir / f"{fasta_path.stem}_proteins.faa"
            with open(final_proteins, "w") as outf:
                for faa in split_dir.glob("*.faa"):
                    with open(faa) as inf:
                        outf.write(inf.read())
            
            # Clean up split directory
            shutil.rmtree(split_dir, ignore_errors=True)
            
            if not final_proteins.exists() or final_proteins.stat().st_size == 0:
                return None
            
            # Compress and delete original
            gz = FileValidator.ensure_gz(final_proteins)
            if gz:
                return gz
            else:
                logger.error(f"Metagenome protein compression failed for {final_proteins}")
                return None
                
        except Exception as e:
            logger.error(f"Metagenome processing failed for {fasta_path}: {e}")
            # Clean up on error
            if split_dir.exists():
                shutil.rmtree(split_dir, ignore_errors=True)
            return None