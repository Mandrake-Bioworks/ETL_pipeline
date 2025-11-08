#!/usr/bin/env python3
"""ENA Extractor - Using proven search API + browser FASTA download"""
import logging
import requests
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

logger = logging.getLogger(__name__)


class ENAExtractor:
    SEARCH_API = "https://www.ebi.ac.uk/ena/portal/api/search"
    BROWSER_FASTA = "https://www.ebi.ac.uk/ena/browser/api/fasta/{}?download=true"

    def __init__(self, config, db):
        self.config = config
        self.db = db

        self.base_dir = Path(config["paths"]["base_data"]) / "ena"
        self.base_dir.mkdir(parents=True, exist_ok=True)

        self.timeout = int(config["processing"]["download_timeout"])
        self.max_retries = int(config["processing"]["max_retries"])
        self.batch_size = int(config["sources"]["ena"]["batch_size"])

        # Kingdoms to fetch
        self.kingdoms = config["sources"]["ena"].get("kingdoms", ["bacteria", "archaea", "viral"])
        
        self.cursor_path = self.base_dir / ".ena_cursor"
        self.metadata = {}

        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "Mandrake-ETL-ENA/2.0",
            "Accept": "application/json",
        })

    def download_batch(self, batch_size, seen_ids):
        """Download a batch of assemblies across all kingdoms"""
        candidates = []
        
        for kingdom in self.kingdoms:
            if len(candidates) >= batch_size:
                break
            
            taxon_id = self._get_taxon_id(kingdom)
            if not taxon_id:
                continue
                
            # Search for complete genomes/chromosomes with more specific query
            query = f'tax_tree({taxon_id}) AND (assembly_level="complete genome" OR assembly_level="chromosome")'
            
            logger.info(f"ENA: Searching {kingdom} with query: {query}")
            assemblies = self._search_assemblies(query, limit=batch_size * 3)
            logger.info(f"ENA: Found {len(assemblies)} assemblies for {kingdom}")
            
            # Filter for complete genomes/chromosomes
            filtered_count = 0
            for assembly in assemblies:
                acc = assembly.get("accession")
                version = assembly.get("version")
                
                if not acc:
                    continue
                
                # Build full accession with version
                if version and not acc.endswith(f".{version}"):
                    full_acc = f"{acc}.{version}"
                else:
                    full_acc = acc
                
                # Filter by assembly level (should already be filtered by query, but double-check)
                assembly_level = (assembly.get("assembly_level") or "").lower()
                if assembly_level not in ["complete genome", "chromosome", "complete"]:
                    continue
                
                filtered_count += 1
                acc_root = acc.split(".", 1)[0]
                
                # Check if already seen (check both with and without version)
                if full_acc in seen_ids or acc in seen_ids or acc_root in seen_ids:
                    logger.debug(f"Skipping {full_acc}: already in seen_ids")
                    continue
                    
                if self.db.entry_exists(full_acc) or self.db.entry_exists(acc) or self.db.entry_exists(acc_root):
                    logger.debug(f"Skipping {full_acc}: already in database")
                    continue
                
                # Store metadata with full accession
                species = self._extract_species(assembly.get("scientific_name"))
                self.metadata[full_acc] = {
                    "kingdom": kingdom,
                    "species": species
                }
                self.metadata[acc] = self.metadata[full_acc]
                self.metadata[acc_root] = self.metadata[full_acc]
                
                # Update assembly dict to use full accession
                assembly["accession"] = full_acc
                
                logger.debug(f"Selected {full_acc}: {species or 'unknown species'} ({assembly_level})")
                candidates.append(assembly)
                
                if len(candidates) >= batch_size:
                    break
            
            logger.info(f"ENA: {kingdom} - {filtered_count} complete genomes found, {len(candidates)} candidates selected")
        
        logger.info(f"ENA: Found {len(candidates)} new assemblies to download")
        
        if not candidates:
            logger.warning("ENA: No candidates found. This might indicate an API issue.")
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
            logger.warning(f"ENA: Failed to download {len(failed)} assemblies: {', '.join(failed[:10])}")
        return downloaded

    def get_metadata(self, accession: str):
        """Get metadata for an accession"""
        return self.metadata.get(accession, {})

    def _get_taxon_id(self, kingdom: str) -> str | None:
        """Map kingdom to NCBI taxonomy ID"""
        mapping = {
            "bacteria": "2",
            "archaea": "2157",
            "viral": "10239"
        }
        return mapping.get(kingdom.lower())

    def _search_assemblies(self, query: str, limit: int = 100):
        """Search ENA for assemblies with better error handling"""
        all_results = []
        offset = 0
        page_size = min(1000, limit)
        
        while len(all_results) < limit:
            params = {
                "result": "assembly",
                "query": query,
                "fields": "accession,scientific_name,tax_division,assembly_level,version",
                "limit": page_size,
                "offset": offset,
                "format": "json"
            }
            
            try:
                response = self.session.get(
                    self.SEARCH_API, 
                    params=params, 
                    timeout=self.timeout
                )
                
                if response.status_code == 200:
                    results = response.json() or []
                    if not results:
                        break
                    
                    all_results.extend(results)
                    
                    if len(results) < page_size:
                        break
                    
                    offset += page_size
                else:
                    logger.error(f"ENA search failed: HTTP {response.status_code}")
                    logger.error(f"Response: {response.text[:500]}")
                    break
                
            except Exception as e:
                logger.error(f"ENA search error: {e}")
                break
        
        return all_results[:limit]

    def _download_assembly(self, assembly: dict):
        """Download assembly using ENA sequence retrieval API with multiple fallback endpoints"""
        acc = assembly.get("accession")
        if not acc:
            return None
        
        # Use ENA's sequence retrieval API with multiple endpoint options
        acc_base = acc.split(".", 1)[0]
        
        # Try different endpoints in order of reliability
        endpoints = [
            # Primary endpoint - direct FASTA download
            f"https://www.ebi.ac.uk/ena/browser/api/fasta/{acc_base}?download=true",
            # Alternative endpoint 1 - ENA data view with FASTA format
            f"https://www.ebi.ac.uk/ena/data/view/{acc_base}&display=fasta&download=fasta",
            # Alternative endpoint 2 - ENA data view without download flag
            f"https://www.ebi.ac.uk/ena/data/view/{acc_base}&display=fasta"
        ]
        
        dest = self.base_dir / f"{acc}.fasta"
        
        # Check if already exists
        if dest.exists() and dest.stat().st_size > 1000:
            logger.info(f"Already exists: {acc}")
            return dest
        
        for endpoint in endpoints:
            for attempt in range(self.max_retries):
                try:
                    logger.info(f"Trying endpoint {endpoints.index(endpoint) + 1}/{len(endpoints)} for {acc}: {endpoint}")
                    
                    # Use proper headers for FASTA download
                    headers = {
                        "User-Agent": "Mandrake-ETL-ENA/2.0",
                        "Accept": "text/x-fasta,text/plain,*/*",
                    }
                    
                    response = self.session.get(endpoint, headers=headers, stream=True, timeout=self.timeout)
                    
                    if response.status_code == 200:
                        # Download to temporary file first
                        temp_dest = dest.with_suffix(".part")
                        bytes_written = 0
                        
                        with open(temp_dest, "wb") as f:
                            for chunk in response.iter_content(chunk_size=1024 * 1024):
                                if chunk:
                                    f.write(chunk)
                                    bytes_written += len(chunk)
                        
                        # Validate it's actually FASTA format (starts with >)
                        is_fasta = False
                        try:
                            with open(temp_dest, "rb") as f:
                                first_bytes = f.read(10)
                                if first_bytes.startswith(b">"):
                                    is_fasta = True
                                else:
                                    logger.warning(f"Invalid FASTA format from {endpoint} - starts with: {first_bytes}")
                        except Exception as e:
                            logger.warning(f"Error validating FASTA format: {e}")
                        
                        if is_fasta and temp_dest.stat().st_size > 1000:
                            # Valid FASTA file - move to final location
                            temp_dest.rename(dest)
                            logger.info(f"✓ Downloaded: {acc} ({bytes_written:,} bytes) from endpoint {endpoints.index(endpoint) + 1}")
                            return dest
                        else:
                            # Invalid content or too small
                            logger.warning(f"Downloaded content invalid or too small for {acc}: {temp_dest.stat().st_size} bytes")
                            if temp_dest.exists():
                                temp_dest.unlink()
                            
                            if attempt < self.max_retries - 1:
                                time.sleep(2)
                                continue
                    else:
                        logger.warning(f"HTTP {response.status_code} from {endpoint}")
                        if attempt < self.max_retries - 1:
                            time.sleep(2)
                            continue
                
                except requests.exceptions.Timeout:
                    logger.warning(f"Timeout downloading {acc} from {endpoint} (attempt {attempt + 1})")
                    if attempt < self.max_retries - 1:
                        time.sleep(2)
                except Exception as e:
                    logger.warning(f"Download attempt {attempt + 1} failed for {acc} from {endpoint}: {e}")
                    if attempt < self.max_retries - 1:
                        time.sleep(2)
            
            # If we successfully downloaded from this endpoint, break the endpoint loop
            if dest.exists() and dest.stat().st_size > 1000:
                break
        
        if not dest.exists() or dest.stat().st_size <= 1000:
            logger.error(f"Failed to download {acc} after trying {len(endpoints)} endpoints with {self.max_retries} attempts each")
            return None
        
        return dest

    @staticmethod
    def _extract_species(scientific_name: str | None):
        """Extract species from scientific name"""
        if not scientific_name:
            return None
        
        tokens = scientific_name.strip().split()
        if len(tokens) >= 2:
            return f"{tokens[0]} {tokens[1]}"
        return tokens[0] if tokens else None