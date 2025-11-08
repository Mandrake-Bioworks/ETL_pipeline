#!/usr/bin/env python3
"""NCBI Extractor - Bacteria, Archaea, Viruses"""
import re
import requests
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

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
        self.metadata = {}  # accession or root -> {'kingdom': ..., 'species': ...}

    @staticmethod
    def _root(acc: str) -> str:
        # "GCF_016406305.1" -> "GCF_016406305"
        return acc.split('.', 1)[0] if acc else acc

    def download_batch(self, batch_size, seen_ids):
        """Download a batch of genomes."""
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
                            acc_full = file_path.name.split('_', 1)[0]  # e.g., GCF_016406305.1
                            seen_ids.add(acc_full)
                    except Exception:
                        pass

        logger.info(f"NCBI: Downloaded {len(downloaded)} genomes")
        return downloaded

    def _get_genomes(self, kingdom, limit, seen_ids):
        """Get genome list from assembly summary."""
        url = f"https://ftp.ncbi.nlm.nih.gov/genomes/refseq/{kingdom}/assembly_summary.txt"
        summary_path = self.base_dir / f"assembly_summary_{kingdom}.txt"

        try:
            r = requests.get(url, timeout=self.timeout)
            summary_path.write_text(r.text)
        except Exception:
            if not summary_path.exists():
                return []

        genomes = []
        with open(summary_path) as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue

                parts = line.split('\t')
                if len(parts) < 20:
                    continue

                acc_full = parts[0]           # assembly_accession (with version)
                acc_root = self._root(acc_full)
                organism_name = parts[7]      # organism_name
                level = parts[11]             # assembly_level
                ftp_path = parts[19]          # ftp_path

                if acc_full in seen_ids or acc_root in seen_ids or self.db.entry_exists(acc_root) or self.db.entry_exists(acc_full):
                    continue

                if level in ['Complete Genome', 'Chromosome'] and ftp_path != 'na':
                    genomes.append((acc_full, ftp_path, kingdom))
                    # Normalize organism_name → "Genus species"
                    sp = None
                    if organism_name:
                        toks = organism_name.strip().split()
                        if len(toks) >= 2:
                            sp = f"{toks[0]} {toks[1]}"
                        elif toks:
                            sp = toks[0]
                    meta = {'kingdom': kingdom, 'species': sp}
                    # store under BOTH keys
                    self.metadata[acc_full] = meta
                    self.metadata[acc_root] = meta

                if len(genomes) >= limit:
                    break

        return genomes

    def _download_genome(self, genome_info):
        """Download single genome."""
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

            # Ensure metadata stays under BOTH keys
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
