#!/usr/bin/env python3
"""
FULL PARALLEL ETL with immediate cleanup after S3 upload
FIXED: Continues through high-duplication regions instead of stopping early
"""
import sys
import yaml
import logging
import logging.handlers
import re
from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock
from Bio import SeqIO

sys.path.insert(0, 'src')

from src.utils.database import Database
from src.utils.helpers import FileValidator, S3Manager, ProdigalRunner, DiskManager
from src.extractors.ncbi_extractor import NCBIExtractor
from src.extractors.ena_extractor import ENAExtractor
from src.extractors.mgnify_extractor import MGnifyExtractor


class ETLPipeline:
    def __init__(self, config_path='etl_config.yaml'):
        with open(config_path) as f:
            self.config = yaml.safe_load(f)

        # logging
        logs_dir = Path(self.config['paths']['logs'])
        logs_dir.mkdir(parents=True, exist_ok=True)
        log_file = logs_dir / f'pipeline_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'
        root = logging.getLogger()
        root.setLevel(logging.INFO)
        sh = logging.StreamHandler()
        sh.setLevel(logging.INFO)
        root.addHandler(sh)
        try:
            fh = logging.handlers.RotatingFileHandler(log_file, maxBytes=10*1024*1024, backupCount=3)
            fh.setLevel(logging.INFO)
            root.addHandler(fh)
        except OSError:
            root.warning("File logging disabled due to filesystem error")

        self.db = Database(self.config)
        self.s3 = S3Manager(self.config)
        self.validator = FileValidator()
        self.prodigal = ProdigalRunner()
        self.disk = DiskManager()

        self.seen_ids = set()
        self.seen_hashes = self.db.get_existing_hashes()

        self.hash_lock = Lock()
        self.log_lock = Lock()

        self.num_workers = int(self.config.get('processing', {}).get('workers', 4))
        self.min_free_gb = int(self.config.get("processing", {}).get("min_free_gb", 10))

        logging.info(f"Pipeline initialized ({len(self.seen_hashes)} existing hashes) | workers={self.num_workers}")

    def run(self):
        logging.info("="*70)
        logging.info("STARTING ETL PIPELINE")
        logging.info("="*70)

        extractors = {
            'ncbi': NCBIExtractor(self.config, self.db),
            'ena' : ENAExtractor(self.config, self.db),
            'mgnify': MGnifyExtractor(self.config, self.db)
        }

        for source in self.config['sources']['order']:
            if source not in extractors or not self.config['sources'].get(source, {}).get('enabled', True):
                continue
            logging.info(f"\n{'='*70}")
            logging.info(f"Processing: {source.upper()}")
            logging.info(f"{'='*70}")
            self._process_source(source, extractors[source])

        self._print_stats()
        self.db.close()
        logging.info("\n✅ Pipeline Complete!")

    def _process_source(self, source, extractor):
        """
        Process a source with proper handling of high-duplication regions.
        FIXED: No longer exits on first empty batch - continues through duplicates.
        """
        batch_size = self.config['sources'][source]['batch_size']
        limit = self.config['sources'][source]['limit']
        total_processed = 0
        batch_num = 0
        consecutive_empty = 0  # Track consecutive empty batches
        is_metagenome = (source == 'mgnify')

        while total_processed < limit:
            batch_num += 1
            remaining = limit - total_processed
            current_batch = min(batch_size, remaining)

            # Disk guard at batch start
            purge_roots = [Path(self.config['paths']['temp']), Path(self.config['paths']['base_data'])]
            if not self.disk.ensure_free_space(self.min_free_gb, purge_paths=purge_roots):
                logging.error("Low disk: cannot ensure %d GB free. Aborting batch.", self.min_free_gb)
                return

            logging.info(f"\n--- Batch {batch_num}: Extracting {current_batch}... ---")
            files = extractor.download_batch(current_batch, self.seen_ids)
            
            # FIXED: Don't break immediately on empty batch
            if not files:
                # Check if extractor reports all sources exhausted
                if hasattr(extractor, 'search_exhausted'):
                    if all(extractor.search_exhausted.values()):
                        logging.info(f"{source.upper()}: All search sources exhausted")
                        break
                
                # Allow multiple consecutive empty batches before giving up
                # (High duplication can cause empty batches with unique data at higher offsets)
                consecutive_empty += 1
                max_empty = 5  # Allow 5 consecutive empty batches
                
                if consecutive_empty >= max_empty:
                    logging.info(f"{source.upper()}: {max_empty} consecutive empty batches, stopping")
                    break
                
                logging.info(f"{source.upper()}: Empty batch {consecutive_empty}/{max_empty}, trying next offset...")
                continue  # Try next batch at higher offset
            
            # Reset consecutive empty counter on successful batch
            consecutive_empty = 0

            logging.info(f"Extracted: {len(files)} files")
            processed_files = []
            with ThreadPoolExecutor(max_workers=self.num_workers) as ex:
                futs = [ex.submit(self._etl_single_file, source, p, is_metagenome, extractor) for p in files]
                for fut in as_completed(futs):
                    ok, p = fut.result()
                    if ok:
                        processed_files.append(p)
                        total_processed += 1

            logging.info(f"Processed: {len(processed_files)} files")
            # Note: files are now cleaned immediately after upload, but we still clean up any remnants
            self.disk.cleanup_files(processed_files)
            
            # REMOVED: Early exit on partial batch
            # This was causing premature termination in high-duplication scenarios
            # if len(files) < current_batch:
            #     break

        logging.info(f"\n{source.upper()} Total: {total_processed} sequences processed")

    @staticmethod
    def _acc_forms_from_filename(p: Path):
        name = p.name
        m = re.search(r'(GC[AF]_\d+(?:\.\d+)?)', name)
        if m:
            acc_full = m.group(1)
            acc_root = acc_full.split('.', 1)[0]
            return acc_full, acc_root
        stem = p.stem
        if re.fullmatch(r'GC[AF]_\d+(?:\.\d+)?', stem):
            return stem, stem.split('.', 1)[0]
        m2 = re.match(r'(GC[AF]_\d+)', stem)
        if m2:
            token = m2.group(1)
            return token, token
        acc_full = stem
        acc_root = stem.split('.', 1)[0]
        return acc_full, acc_root

    def _etl_single_file(self, source, file_path, is_metagenome, extractor):
        try:
            # Disk guard per file
            purge_roots = [Path(self.config['paths']['temp'])]
            if not self.disk.ensure_free_space(self.min_free_gb, purge_paths=purge_roots):
                with self.log_lock:
                    logging.error("Low disk during ETL of %s — skipping.", file_path.name)
                return False, file_path

            # 1) VALIDATE
            valid_path = self.validator.validate_and_fix(file_path)
            if not valid_path:
                with self.log_lock:
                    logging.warning(f"Validation failed: {file_path.name}")
                return False, file_path

            # 1a) METADATA lookup
            acc_full, acc_root = self._acc_forms_from_filename(file_path)
            meta = {}
            if hasattr(extractor, 'get_metadata'):
                meta = (extractor.get_metadata(acc_full) or extractor.get_metadata(acc_root) or {})
            kingdom = meta.get('kingdom') if not is_metagenome else None
            origin  = meta.get('origin')  if is_metagenome  else None
            species_meta = meta.get('species')

            # 1b) species from header
            species_parsed = None
            name_lc = file_path.name.lower()
            if "cds_from_genomic" not in name_lc and "cds-" not in name_lc:
                species_parsed, _ = self.validator.parse_species_family(valid_path)

            def _clean_species(s):
                if not s:
                    return None
                s = s.strip()
                if "[" in s or "]" in s or len(s.split()) < 2:
                    return None
                return s

            species = _clean_species(species_parsed) or _clean_species(species_meta)

            # 2) DEDUP by sequence hash
            seq_hash = self.validator.compute_sequence_hash(valid_path)
            if not seq_hash:
                with self.log_lock:
                    logging.warning(f"Hash failed: {file_path.name}")
                return False, file_path

            with self.hash_lock:
                if seq_hash in self.seen_hashes:
                    logging.info(f"Duplicate (in-mem): {acc_full or file_path.name}")
                    return False, file_path
                self.seen_hashes.add(seq_hash)

            # 3) FILTER (with stats collection)
            if is_metagenome:
                filtered_path, filter_stats = self._filter_metagenome(valid_path)
                # Record filtering stats
                if filter_stats:
                    self.db.insert_filtering_stats(
                        source=source,
                        accession=(acc_full or acc_root or file_path.stem),
                        total_contigs=filter_stats['total'],
                        contigs_kept=filter_stats['kept'],
                        contigs_removed=filter_stats['removed']
                    )
            else:
                filtered_path = valid_path

            # 4) CALCULATE TOTAL BP BEFORE COMPRESSION
            total_bp = sum(len(r.seq) for r in SeqIO.parse(filtered_path, 'fasta'))

            # 5) CONVERT
            if is_metagenome:
                work_dir = Path(self.config['paths']['temp']) / (acc_full or acc_root or file_path.stem)
                protein_path_gz = self.prodigal.process_metagenome(filtered_path, work_dir)
            else:
                protein_path_gz = self.prodigal.predict_proteins(filtered_path, is_metagenome=False)

            if not protein_path_gz:
                with self.log_lock:
                    logging.warning(f"Protein prediction failed: {acc_full or file_path.name}")
                return False, file_path

            genome_upload_path_gz = self.validator.ensure_gz(filtered_path) or filtered_path

            # 6) LOAD (S3 + DB) with immediate cleanup
            ok = self._load_and_cleanup(
                source=source,
                accession=(acc_full or acc_root or file_path.stem),
                genome_path_gz=genome_upload_path_gz,
                protein_path_gz=protein_path_gz,
                seq_hash=seq_hash,
                species=species,
                kingdom=kingdom,
                origin=origin,
                total_bp=total_bp  # Pass pre-calculated value instead of file path
            )
            
            # cleanup per-analysis work dir for metagenomes
            if ok and is_metagenome:
                work_dir = Path(self.config['paths']['temp']) / (acc_full or acc_root or file_path.stem)
                self.disk.cleanup_directory(work_dir)
            
            return ok, file_path

        except Exception as e:
            with self.log_lock:
                logging.error(f"ETL error for {file_path.name}: {e}")
            return False, file_path

    def _filter_metagenome(self, file_path):
        """Filter metagenome and return stats"""
        min_length = int(self.config['filtering']['metagenomes']['min_contig_length'])
        filtered = file_path.with_suffix('.filtered.fasta')
        
        total_contigs = 0
        kept_contigs = 0
        kept_records = []
        
        for r in SeqIO.parse(file_path, 'fasta'):
            total_contigs += 1
            if len(r.seq) > min_length:
                kept_records.append(r)
                kept_contigs += 1
        
        removed_contigs = total_contigs - kept_contigs
        
        stats = {
            'total': total_contigs,
            'kept': kept_contigs,
            'removed': removed_contigs
        }
        
        if kept_records:
            from Bio import SeqIO as _SeqIO
            _SeqIO.write(kept_records, filtered, 'fasta')
            return filtered, stats
        
        return file_path, stats

    def _load_and_cleanup(self, source, accession, genome_path_gz, protein_path_gz, 
                         seq_hash, species, kingdom, origin, total_bp):
        """Load to S3 and DB, with immediate cleanup after each upload"""
        try:
            # Upload genome and cleanup immediately
            genome_s3 = self.s3.upload_genome(genome_path_gz, source, accession)
            if genome_s3:
                # Delete genome file immediately after successful upload
                try:
                    if genome_path_gz.exists():
                        genome_path_gz.unlink()
                        logging.info(f"Cleaned up genome file: {genome_path_gz.name}")
                except Exception as e:
                    logging.warning(f"Failed to cleanup genome {genome_path_gz}: {e}")
            else:
                logging.error(f"Genome upload failed: {genome_path_gz}")
                return False

            # Upload proteins and cleanup immediately
            protein_s3 = self.s3.upload_proteins(protein_path_gz, source, accession)
            if protein_s3:
                # Delete protein file immediately after successful upload
                try:
                    if protein_path_gz.exists():
                        protein_path_gz.unlink()
                        logging.info(f"Cleaned up protein file: {protein_path_gz.name}")
                except Exception as e:
                    logging.warning(f"Failed to cleanup protein {protein_path_gz}: {e}")
            else:
                logging.error(f"Protein upload failed: {protein_path_gz}")
                return False

            # Insert into database - total_bp is now pre-calculated
            inserted, reason = self.db.insert_entry(
                source=source,
                accession=accession,
                s3_genome_path=genome_s3,
                s3_protein_path=protein_s3,
                sequence_hash=seq_hash,
                total_bp=total_bp,
                species=species,
                kingdom=kingdom,
                origin=origin
            )

            if not inserted:
                with self.log_lock:
                    if reason == "hash_conflict":
                        logging.info(f"Duplicate (DB race, sequence_hash): {accession}")
                    elif reason == "accession_conflict":
                        logging.info(f"Duplicate (DB race, source+accession): {source}:{accession}")
                    else:
                        logging.info(f"Duplicate (DB race): {source}:{accession}")
                return False

            with self.log_lock:
                logging.info(
                    f"✓ LOADED: {accession} ({total_bp:,} bp)  "
                    f"{'species='+species+' ' if species else ''}"
                    f"{'kingdom='+kingdom+' ' if kingdom else ''}"
                    f"{'origin='+origin if origin else ''}"
                )
            return True

        except Exception as e:
            with self.log_lock:
                logging.error(f"Load failed for {accession}: {e}", exc_info=True)
            return False

    def _print_stats(self):
        logging.info("\n" + "="*70)
        logging.info("FINAL STATISTICS")
        logging.info("="*70)
        stats = self.db.get_stats()
        for row in stats:
            source, entries, bp, species = row
            logging.info(f"{source:15s}: {entries:6d} entries, {bp:12,d} bp, {species:4d} species")


def main():
    pipeline = ETLPipeline()
    pipeline.run()


if __name__ == '__main__':
    main()
