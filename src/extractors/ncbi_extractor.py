#!/usr/bin/env python3
"""NCBI Extractor - Enhanced with genome_rep filtering for full genomes only"""
import re
import requests
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
from datetime import datetime

logger = logging.getLogger(__name__)


class NCBIExtractor:
    def __init__(self, config, db):
        self.config = config
        self.db = db
        self.base_dir = Path(config['paths']['base_data']) / 'ncbi'
        self.base_dir.mkdir(parents=True, exist_ok=True)
        self.timeout = config['processing']['download_timeout']
        self.kingdoms = config['sources']['ncbi'].get('kingdoms', ['bacteria', 'archaea'])
        self.batch_size = config['sources']['ncbi']['batch_size']
        self.limit = config['sources']['ncbi']['limit']
        self.metadata = {}
        
        self.force_update_summaries = config['sources']['ncbi'].get('force_update_summaries', False)
        
        # Filter for full genome representation only (not partial)
        self.require_full_genome = config['sources']['ncbi'].get('require_full_genome', True)
        
        self.existing_accessions = self._load_existing_accessions()
        logger.info(f"NCBI: Found {len(self.existing_accessions)} existing accessions in database")

    def _load_existing_accessions(self):
        """Load all existing NCBI accessions from database"""
        conn = self.db._getconn()
        try:
            with conn.cursor() as cur:
                cur.execute("""
                    SELECT DISTINCT accession 
                    FROM entries 
                    WHERE source = 'ncbi'
                """)
                accessions = set()
                for row in cur.fetchall():
                    acc = row[0]
                    accessions.add(acc)
                    accessions.add(self._root(acc))
                    
                return accessions
        finally:
            self.db._putconn(conn)

    @staticmethod
    def _root(acc: str) -> str:
        """Extract root accession without version"""
        return acc.split('.', 1)[0] if acc else acc
    
    @staticmethod
    def _extract_accession_from_filename(filename: str) -> str:
        """
        Extract full accession from NCBI filename
        Example: "GCF_002271195.1_ASM227119v1_genomic.fna.gz" -> "GCF_002271195.1"
        """
        # Pattern: GCF_XXXXXXXXX.X or GCA_XXXXXXXXX.X
        match = re.match(r'(GC[FA]_\d+\.\d+)', filename)
        if match:
            return match.group(1)
        
        # Fallback: Try to extract from the first underscore-separated segment
        parts = filename.split('_')
        if len(parts) >= 2:
            potential_acc = f"{parts[0]}_{parts[1]}"
            if re.match(r'GC[FA]_\d+\.\d+', potential_acc):
                return potential_acc
        
        return filename.split('_genomic')[0] if '_genomic' in filename else filename

    def mark_as_processed(self, accession: str):
        """Mark an accession as processed"""
        if accession:
            self.existing_accessions.add(accession)
            self.existing_accessions.add(self._root(accession))
            logger.debug(f"Marked as processed: {accession}")

    def download_batch(self, batch_size, seen_ids):
        """Download a batch of genomes"""
        downloaded = []

        for kingdom in self.kingdoms:
            if len(downloaded) >= batch_size:
                break

            genomes = self._get_genomes(kingdom, batch_size - len(downloaded), seen_ids)

            with ThreadPoolExecutor(max_workers=8) as executor:
                futures = {executor.submit(self._download_genome, g): g for g in genomes}

                for future in as_completed(futures):
                    try:
                        file_path = future.result()
                        if file_path:
                            downloaded.append(file_path)
                            
                            # Extract accession properly from filename
                            acc_full = self._extract_accession_from_filename(file_path.name)
                            
                            seen_ids.add(acc_full)
                            seen_ids.add(self._root(acc_full))
                            
                            # Update cache immediately
                            self.mark_as_processed(acc_full)
                            
                            logger.debug(f"Downloaded and cached: {acc_full}")
                    except Exception as e:
                        logger.debug(f"Download exception: {e}")

        logger.info(f"NCBI: Downloaded {len(downloaded)} genomes")
        logger.info(f"NCBI: Cache now has {len(self.existing_accessions)} accessions")
        return downloaded

    def _get_genomes(self, kingdom, limit, seen_ids):
        """Get genome list with enhanced filtering"""
        url = f"https://ftp.ncbi.nlm.nih.gov/genomes/refseq/{kingdom}/assembly_summary.txt"
        summary_path = self.base_dir / f"assembly_summary_{kingdom}.txt"

        # Check if assembly summary exists
        if summary_path.exists() and summary_path.stat().st_size > 1000:
            file_size_mb = summary_path.stat().st_size / (1024 * 1024)
            age_days = (time.time() - summary_path.stat().st_mtime) / 86400
            
            if self.force_update_summaries:
                logger.info(f"NCBI: Force update enabled, re-downloading for {kingdom}...")
                self._download_assembly_summary(url, summary_path)
            else:
                logger.info(f"NCBI: Using existing assembly summary for {kingdom}")
                logger.info(f"      File: {summary_path.name} ({file_size_mb:.1f} MB, {age_days:.0f} days old)")
                file_age = datetime.fromtimestamp(summary_path.stat().st_mtime)
                logger.info(f"      Last modified: {file_age.strftime('%Y-%m-%d %H:%M:%S')}")
        else:
            logger.info(f"NCBI: Assembly summary not found, downloading for {kingdom}...")
            self._download_assembly_summary(url, summary_path)

        if not summary_path.exists() or summary_path.stat().st_size < 1000:
            logger.error(f"NCBI: Assembly summary invalid: {summary_path}")
            return []

        # Parse assembly summary
        logger.info(f"NCBI: Parsing assembly summary for {kingdom}...")
        genomes = []
        skipped_in_cache = 0
        skipped_seen = 0
        skipped_level = 0
        skipped_partial = 0
        total_lines = 0
        
        start_parse = time.time()
        
        with open(summary_path) as f:
            for line in f:
                total_lines += 1
                
                if total_lines % 50000 == 0:
                    logger.info(f"NCBI: Parsed {total_lines:,} lines, found {len(genomes)} candidates")
                
                if line.startswith('#') or not line.strip():
                    continue

                parts = line.split('\t')
                if len(parts) < 20:
                    continue

                # Column indices (0-based)
                acc_full = parts[0]          # Column 1: assembly_accession
                organism_name = parts[7]     # Column 8: organism_name
                level = parts[11]            # Column 12: assembly_level
                genome_rep = parts[13]       # Column 14: genome_rep
                ftp_path = parts[19]         # Column 20: ftp_path
                
                acc_root = self._root(acc_full)

                # Check cache (includes DB + session updates)
                if acc_full in self.existing_accessions or acc_root in self.existing_accessions:
                    skipped_in_cache += 1
                    continue

                # Check seen_ids (this batch)
                if acc_full in seen_ids or acc_root in seen_ids:
                    skipped_seen += 1
                    continue

                # Assembly level filter
                if level not in ['Complete Genome', 'Chromosome'] or ftp_path == 'na':
                    skipped_level += 1
                    continue
                
                # NEW: Genome representation filter
                if self.require_full_genome and genome_rep != 'Full':
                    skipped_partial += 1
                    continue

                # Accept this genome
                genomes.append((acc_full, ftp_path, kingdom))
                
                # Store metadata
                sp = None
                if organism_name:
                    toks = organism_name.strip().split()
                    if len(toks) >= 2:
                        sp = f"{toks[0]} {toks[1]}"
                    elif toks:
                        sp = toks[0]
                meta = {'kingdom': kingdom, 'species': sp}
                self.metadata[acc_full] = meta
                self.metadata[acc_root] = meta

                if len(genomes) >= limit:
                    break

        parse_time = time.time() - start_parse
        logger.info(f"NCBI: Parsed {total_lines:,} lines in {parse_time:.1f}s")
        logger.info(f"NCBI {kingdom}: {len(genomes)} new genomes, "
                   f"skipped {skipped_in_cache} in cache, "
                   f"{skipped_seen} in session, {skipped_level} wrong level, "
                   f"{skipped_partial} partial genomes")
        
        return genomes

    def _download_assembly_summary(self, url: str, dest_path: Path):
        """Download assembly summary file"""
        try:
            start_time = time.time()
            r = requests.get(url, timeout=self.timeout, stream=True)
            total_size = int(r.headers.get('content-length', 0))
            logger.info(f"NCBI: Downloading {total_size / 1_000_000:.1f} MB...")
            
            with open(dest_path, 'wb') as f:
                downloaded = 0
                last_log = 0
                for chunk in r.iter_content(chunk_size=1024*1024):
                    if chunk:
                        f.write(chunk)
                        downloaded += len(chunk)
                        if downloaded - last_log >= 100*1024*1024:
                            logger.info(f"NCBI: Downloaded {downloaded / 1_000_000:.1f} MB...")
                            last_log = downloaded
            
            download_time = time.time() - start_time
            logger.info(f"NCBI: Downloaded ({download_time:.1f}s, "
                       f"{(total_size / download_time / 1_000_000):.1f} MB/s)")
            
        except Exception as e:
            logger.error(f"NCBI: Failed to download: {e}")
            raise

    def _download_genome(self, genome_info):
        """Download single genome"""
        acc_full, ftp_path, kingdom = genome_info
        acc_root = self._root(acc_full)
        https_url = ftp_path.replace('ftp://', 'https://')

        try:
            html = requests.get(https_url + '/', timeout=self.timeout).text
            match = re.search(rf"({re.escape(acc_full)}[^\"]*genomic\.fna\.gz)", html)
            if not match:
                return None

            filename = match.group(1)
            dest = self.base_dir / filename

            if not dest.exists():
                r = requests.get(f"{https_url}/{filename}", stream=True, timeout=self.timeout)
                if r.ok:
                    with open(dest, 'wb') as f:
                        for chunk in r.iter_content(chunk_size=1024*1024):
                            if chunk:
                                f.write(chunk)

            cur = self.metadata.get(acc_full, {'kingdom': kingdom, 'species': None})
            cur['kingdom'] = cur.get('kingdom') or kingdom
            self.metadata[acc_full] = cur
            self.metadata[acc_root] = cur
            return dest if dest.exists() else None

        except Exception as e:
            logger.debug(f"Download failed for {acc_full}: {e}")
            return None

    def get_metadata(self, accession: str):
        return self.metadata.get(accession, {})
