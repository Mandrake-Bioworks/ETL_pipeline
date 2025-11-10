#!/usr/bin/env python3
"""ENA Extractor v6 - Hybrid TSV + Optional Metadata Enrichment"""
import logging
import requests
import pandas as pd
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
import json
import re

logger = logging.getLogger(__name__)


class ENAExtractor:
    TSV_API = "https://www.ebi.ac.uk/ena/browser/api/tsv/textsearch"
    PORTAL_API = "https://www.ebi.ac.uk/ena/portal/api/search"
    BROWSER_FASTA = "https://www.ebi.ac.uk/ena/browser/api/fasta/{}?download=true"
    BROWSER_XML = "https://www.ebi.ac.uk/ena/browser/api/xml/{}"

    def __init__(self, config, db):
        self.config = config
        self.db = db

        self.base_dir = Path(config["paths"]["base_data"]) / "ena"
        self.base_dir.mkdir(parents=True, exist_ok=True)

        self.timeout = int(config["processing"]["download_timeout"])
        self.max_retries = int(config["processing"]["max_retries"])
        self.batch_size = int(config["sources"]["ena"]["batch_size"])

        self.kingdoms = config["sources"]["ena"].get("kingdoms", ["bacteria", "archaea", "viral"])
        
        # Metadata cache files
        self.tsv_cache = self.base_dir / "ena_catalog.tsv"
        self.metadata_cache = self.base_dir / "ena_metadata.json"
        self.metadata = {}

        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "Mandrake-ETL-ENA/6.0-Hybrid"
        })
        
        # STATEFUL: Track position in TSV catalog
        self.tsv_position = 0
        self.tsv_data = None
        self.catalog_exhausted = False
        self.portal_api_available = None  # Will test on first use
        
        # Load existing accessions at startup
        self.existing_accessions = self._load_existing_accessions()
        logger.info(f"ENA: Found {len(self.existing_accessions)} existing accessions in database")
        
        # Load or download TSV catalog
        self._ensure_catalog()
        
        # Load metadata cache if exists
        self._load_metadata_cache()

    def _load_existing_accessions(self):
        """Load all existing ENA accessions from database"""
        conn = self.db._getconn()
        try:
            with conn.cursor() as cur:
                cur.execute("""
                    SELECT DISTINCT accession 
                    FROM entries 
                    WHERE source = 'ena'
                """)
                accessions = {row[0] for row in cur.fetchall()}
                roots = {acc.split('.', 1)[0] for acc in accessions}
                return accessions | roots
        finally:
            self.db._putconn(conn)

    def _ensure_catalog(self):
        """Download TSV catalog if not already present"""
        if self.tsv_cache.exists():
            logger.info(f"ENA: Loading existing catalog from {self.tsv_cache}")
            try:
                self.tsv_data = pd.read_csv(self.tsv_cache, sep="\t", low_memory=False)
                logger.info(f"ENA: Loaded {len(self.tsv_data)} assemblies from cached catalog")
                return
            except Exception as e:
                logger.warning(f"ENA: Failed to load cached catalog: {e}, will re-download")
        
        logger.info("ENA: Downloading TSV catalog (~14K prokaryotic assemblies)...")
        self._download_tsv_catalog()

    def _download_tsv_catalog(self):
        """Download the full TSV catalog from ENA"""
        params = {
            "domain": "genome_assembly",
            "query": "prokaryotic whole genome sequences"
        }
        
        # Use simple headers for TSV endpoint (complex Accept header causes HTTP 406)
        tsv_headers = {
            "User-Agent": "Mandrake-ETL-ENA/6.0-Hybrid"
        }
        
        try:
            response = requests.get(
                self.TSV_API,
                params=params,
                headers=tsv_headers,
                timeout=120  # Longer timeout for full catalog
            )
            
            if response.status_code == 200:
                # Save raw TSV
                self.tsv_cache.write_text(response.text, encoding='utf-8')
                logger.info(f"ENA: Downloaded catalog to {self.tsv_cache}")
                
                # Parse into DataFrame
                self.tsv_data = pd.read_csv(self.tsv_cache, sep="\t", low_memory=False)
                logger.info(f"ENA: Parsed {len(self.tsv_data)} assemblies from catalog")
                
                # Extract basic metadata from descriptions
                self._extract_metadata_from_descriptions()
                
            else:
                logger.error(f"ENA: Failed to download TSV catalog: HTTP {response.status_code}")
                self.tsv_data = pd.DataFrame()
                
        except Exception as e:
            logger.error(f"ENA: Error downloading TSV catalog: {e}")
            self.tsv_data = pd.DataFrame()

    def _extract_metadata_from_descriptions(self):
        """Extract basic metadata from description field"""
        if self.tsv_data is None or self.tsv_data.empty:
            return
        
        logger.info("ENA: Extracting metadata from descriptions...")
        
        for idx, row in self.tsv_data.iterrows():
            try:
                acc = row.get('accession')
                if pd.isna(acc):
                    continue
                
                acc = str(acc).strip()
                acc_root = acc.split('.', 1)[0]
                
                # Extract species from description
                description = str(row.get('description', ''))
                species = self._extract_species_from_description(description)
                
                # Infer kingdom from species name
                kingdom = self._infer_kingdom_from_species(species)
                
                # Store basic metadata
                metadata_entry = {
                    "kingdom": kingdom,
                    "species": species,
                    "description": description,
                    "assembly_level": "unknown",  # Will try to enrich later
                    "scientific_name": species or "Unknown"
                }
                
                self.metadata[acc] = metadata_entry
                self.metadata[acc_root] = metadata_entry
                
            except Exception as e:
                logger.debug(f"ENA: Error extracting metadata for row {idx}: {e}")
                continue
        
        logger.info(f"ENA: Extracted metadata for {len(self.metadata)} accessions")
        
        # Try to enrich metadata via Portal API (might be available on some networks)
        self._try_enrich_metadata()
        
        # Save metadata cache
        self._save_metadata_cache()

    def _extract_species_from_description(self, description):
        """Extract species name from description field"""
        if not description or pd.isna(description):
            return None
        
        # Pattern: "assembly for Species name"
        match = re.search(r'assembly for (.+?)(?:\s+strain|\s+isolate|\s+\d|$)', description, re.IGNORECASE)
        if match:
            species_text = match.group(1).strip()
            # Take first two words as species name
            parts = species_text.split()[:2]
            if len(parts) == 2:
                return f"{parts[0]} {parts[1]}"
            elif len(parts) == 1:
                return parts[0]
        
        return None

    def _infer_kingdom_from_species(self, species):
        """Infer kingdom from species name patterns"""
        if not species:
            return 'bacteria'  # Default for prokaryotic dataset
        
        species_lower = species.lower()
        
        # Viral patterns
        if any(word in species_lower for word in ['virus', 'phage', 'viroid']):
            return 'viral'
        
        # Archaeal patterns
        if any(word in species_lower for word in ['methanobacterium', 'halobacterium', 'thermococcus', 
                                                     'pyrococcus', 'sulfolobus', 'methanococcus']):
            return 'archaea'
        
        # Default to bacteria for prokaryotic dataset
        return 'bacteria'

    def _try_enrich_metadata(self):
        """Attempt to enrich metadata via Portal API (may fail on restricted networks)"""
        if self.portal_api_available is False:
            logger.debug("ENA: Portal API known to be unavailable, skipping enrichment")
            return
        
        # Test with a small sample
        sample_accessions = list(self.metadata.keys())[:10]
        if not sample_accessions:
            return
        
        logger.info("ENA: Testing Portal API for metadata enrichment...")
        
        try:
            # Build query for sample accessions
            acc_query = " OR ".join([f'accession="{acc.split(".", 1)[0]}"' for acc in sample_accessions[:3]])
            
            params = {
                "result": "assembly",
                "query": acc_query,
                "fields": "accession,scientific_name,assembly_level,tax_division,base_count",
                "format": "json",
                "dataPortal": "ena"
            }
            
            response = self.session.get(self.PORTAL_API, params=params, timeout=30)
            
            if response.status_code == 200:
                data = response.json()
                if data:
                    logger.info(f"ENA: Portal API available! Will enrich metadata for all assemblies...")
                    self.portal_api_available = True
                    self._enrich_all_metadata()
                else:
                    logger.info("ENA: Portal API returned no data, using description-based metadata")
                    self.portal_api_available = False
            else:
                logger.info(f"ENA: Portal API unavailable (HTTP {response.status_code}), using description-based metadata")
                self.portal_api_available = False
                
        except Exception as e:
            logger.info(f"ENA: Portal API unavailable ({type(e).__name__}), using description-based metadata")
            self.portal_api_available = False

    def _enrich_all_metadata(self):
        """Enrich metadata for all accessions using Portal API in batches"""
        if not self.tsv_data is not None:
            return
        
        # Get unique accessions
        accessions = [str(acc).split('.', 1)[0] for acc in self.metadata.keys() if '.' in str(acc)]
        unique_accessions = list(set(accessions))[:1000]  # Limit to first 1000 to avoid long waits
        
        logger.info(f"ENA: Enriching metadata for up to {len(unique_accessions)} assemblies...")
        
        batch_size = 50
        enriched_count = 0
        
        for i in range(0, len(unique_accessions), batch_size):
            batch = unique_accessions[i:i+batch_size]
            
            try:
                acc_query = " OR ".join([f'accession="{acc}"' for acc in batch])
                
                params = {
                    "result": "assembly",
                    "query": acc_query,
                    "fields": "accession,scientific_name,assembly_level,tax_division,base_count",
                    "format": "json",
                    "dataPortal": "ena"
                }
                
                response = self.session.get(self.PORTAL_API, params=params, timeout=60)
                
                if response.status_code == 200:
                    data = response.json()
                    for item in data:
                        acc = item.get('accession')
                        if acc and acc in self.metadata:
                            # Update with richer metadata
                            self.metadata[acc].update({
                                "assembly_level": item.get('assembly_level', 'unknown'),
                                "scientific_name": item.get('scientific_name', self.metadata[acc]['species']),
                                "tax_division": item.get('tax_division', ''),
                                "kingdom": self._infer_kingdom_from_tax_division(item.get('tax_division', ''))
                            })
                            enriched_count += 1
                
                time.sleep(0.5)  # Rate limiting
                
            except Exception as e:
                logger.debug(f"ENA: Error enriching batch {i//batch_size}: {e}")
                continue
        
        logger.info(f"ENA: Enriched metadata for {enriched_count} assemblies")

    def _infer_kingdom_from_tax_division(self, tax_division):
        """Infer kingdom from taxonomy division"""
        tax_div = str(tax_division).lower()
        
        if 'bac' in tax_div or 'pro' in tax_div:
            return 'bacteria'
        elif 'arch' in tax_div:
            return 'archaea'
        elif 'vir' in tax_div:
            return 'viral'
        else:
            return 'bacteria'

    def _load_metadata_cache(self):
        """Load metadata from JSON cache"""
        if self.metadata_cache.exists():
            try:
                with open(self.metadata_cache, 'r') as f:
                    self.metadata = json.load(f)
                logger.info(f"ENA: Loaded metadata cache for {len(self.metadata)} accessions")
            except Exception as e:
                logger.warning(f"ENA: Failed to load metadata cache: {e}")

    def _save_metadata_cache(self):
        """Save metadata to JSON cache for quick reloads"""
        try:
            with open(self.metadata_cache, 'w') as f:
                json.dump(self.metadata, f, indent=2)
            logger.info(f"ENA: Saved metadata cache to {self.metadata_cache}")
        except Exception as e:
            logger.error(f"ENA: Failed to save metadata cache: {e}")

    def download_batch(self, batch_size, seen_ids):
        """Download a batch of assemblies from TSV catalog"""
        if self.tsv_data is None or self.tsv_data.empty:
            logger.warning("ENA: TSV catalog not available")
            self.catalog_exhausted = True
            return []
        
        if self.catalog_exhausted:
            logger.info("ENA: Catalog exhausted - no more data available")
            return []
        
        candidates = []
        skipped_in_db = 0
        skipped_seen = 0
        skipped_kingdom = 0
        
        # Iterate through TSV starting from current position
        total_rows = len(self.tsv_data)
        
        while len(candidates) < batch_size and self.tsv_position < total_rows:
            row = self.tsv_data.iloc[self.tsv_position]
            self.tsv_position += 1
            
            try:
                # Get accession
                acc = row.get('accession')
                if pd.isna(acc):
                    continue
                
                acc = str(acc).strip()
                acc_root = acc.split('.', 1)[0]
                
                # Filter by kingdom if specified
                if self.kingdoms:
                    kingdom = self.metadata.get(acc, {}).get('kingdom', 'bacteria')
                    if kingdom not in self.kingdoms:
                        skipped_kingdom += 1
                        continue
                
                # Check if already processed
                if (acc in self.existing_accessions or 
                    acc_root in self.existing_accessions):
                    skipped_in_db += 1
                    continue
                
                if acc in seen_ids or acc_root in seen_ids:
                    skipped_seen += 1
                    continue
                    
                if self.db.entry_exists(acc) or self.db.entry_exists(acc_root):
                    skipped_in_db += 1
                    self.existing_accessions.add(acc)
                    self.existing_accessions.add(acc_root)
                    continue
                
                # Add to candidates
                assembly = {
                    "accession": acc,
                    "scientific_name": self.metadata.get(acc, {}).get('scientific_name'),
                    "assembly_level": self.metadata.get(acc, {}).get('assembly_level', 'unknown')
                }
                candidates.append(assembly)
                
            except Exception as e:
                logger.debug(f"ENA: Error processing row {self.tsv_position}: {e}")
                continue
        
        # Check if catalog exhausted
        if self.tsv_position >= total_rows:
            self.catalog_exhausted = True
            logger.info("ENA: Reached end of catalog")
        
        logger.info(f"ENA: Position {self.tsv_position}/{total_rows}, "
                   f"found {len(candidates)} new candidates, "
                   f"skipped {skipped_in_db} in DB, {skipped_seen} in session, "
                   f"{skipped_kingdom} wrong kingdom")
        
        if not candidates:
            logger.warning("ENA: No new candidates found")
            return []
        
        # Download in parallel
        downloaded = []
        failed = []
        with ThreadPoolExecutor(max_workers=8) as ex:
            futures = {ex.submit(self._download_assembly, a): a for a in candidates}
            for fut in as_completed(futures):
                try:
                    file_path = fut.result()
                    if file_path:
                        acc = file_path.stem.split(".")[0]
                        seen_ids.add(acc)
                        seen_ids.add(acc.split(".", 1)[0])
                        downloaded.append(file_path)
                    else:
                        assembly = futures[fut]
                        failed.append(assembly.get("accession", "unknown"))
                except Exception as e:
                    logger.error(f"Download thread error: {e}")

        logger.info(f"ENA: Downloaded {len(downloaded)}/{len(candidates)} assemblies")
        if failed:
            logger.warning(f"ENA: Failed to download {len(failed)} assemblies")
        return downloaded

    def get_metadata(self, accession: str):
        """Get metadata for an accession"""
        return self.metadata.get(accession, {})

    def _download_assembly(self, assembly: dict):
        """Download assembly using ENA browser API"""
        acc = assembly.get("accession")
        if not acc:
            return None
        
        acc_base = acc.split(".", 1)[0]
        
        # Primary endpoint (more reliable)
        endpoints = [
            f"https://www.ebi.ac.uk/ena/browser/api/fasta/{acc_base}?download=true",
            f"https://www.ebi.ac.uk/ena/data/view/{acc_base}&display=fasta&download=fasta",
            f"https://www.ebi.ac.uk/ena/data/view/{acc_base}&display=fasta"
        ]
        
        dest = self.base_dir / f"{acc}.fasta"
        
        if dest.exists() and dest.stat().st_size > 1000:
            logger.debug(f"Already exists: {acc}")
            return dest
        
        for endpoint in endpoints:
            for attempt in range(self.max_retries):
                try:
                    headers = {
                        "User-Agent": "Mandrake-ETL-ENA/6.0-Hybrid",
                        "Accept": "text/x-fasta,text/plain,*/*",
                    }
                    
                    response = self.session.get(endpoint, headers=headers, stream=True, timeout=self.timeout)
                    
                    if response.status_code == 200:
                        temp_dest = dest.with_suffix(".part")
                        bytes_written = 0
                        
                        with open(temp_dest, "wb") as f:
                            for chunk in response.iter_content(chunk_size=1024 * 1024):
                                if chunk:
                                    f.write(chunk)
                                    bytes_written += len(chunk)
                        
                        # Validate FASTA format
                        is_fasta = False
                        try:
                            with open(temp_dest, "rb") as f:
                                first_bytes = f.read(10)
                                if first_bytes.startswith(b">"):
                                    is_fasta = True
                        except Exception:
                            pass
                        
                        if is_fasta and temp_dest.stat().st_size > 1000:
                            temp_dest.rename(dest)
                            logger.info(f"✓ Downloaded: {acc} ({bytes_written:,} bytes)")
                            return dest
                        else:
                            if temp_dest.exists():
                                temp_dest.unlink()
                            
                            if attempt < self.max_retries - 1:
                                time.sleep(2)
                                continue
                    else:
                        if attempt < self.max_retries - 1:
                            time.sleep(2)
                            continue
                
                except requests.exceptions.Timeout:
                    if attempt < self.max_retries - 1:
                        time.sleep(2)
                except Exception as e:
                    if attempt < self.max_retries - 1:
                        time.sleep(2)
            
            if dest.exists() and dest.stat().st_size > 1000:
                break
        
        if not dest.exists() or dest.stat().st_size <= 1000:
            logger.error(f"Failed to download {acc}")
            return None
        
        return dest

    @staticmethod
    def _extract_species(scientific_name):
        """Extract species from scientific name"""
        if not scientific_name or pd.isna(scientific_name):
            return None
        
        scientific_name = str(scientific_name).strip()
        tokens = scientific_name.split()
        
        if len(tokens) >= 2:
            return f"{tokens[0]} {tokens[1]}"
        return tokens[0] if tokens else None
