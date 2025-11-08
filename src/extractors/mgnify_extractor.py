#!/usr/bin/env python3
"""MGnify Extractor - ONLY processed genomic reads (FASTQ.fasta[.gz]) with rate limiting + origin metadata"""
import time
import json
import gzip
import logging
import requests
from pathlib import Path

logger = logging.getLogger(__name__)


class MGnifyExtractor:
    """
    Download ONLY processed genomic reads from MGnify analyses:
      ✅ IDs that match *FASTQ.fasta or *FASTQ.fasta.gz
      ❌ Skip predicted ORFs (.ffn), proteins (.faa), contigs/assembly .fasta

    Config (optional):
      sources.mgnify.environments:          list[str]  (default: ['soil','marine','freshwater','plant','gut','sediment'])
      sources.mgnify.analyses_per_study:    int        (default: 2)
      sources.mgnify.delay_seconds:         float      (default: 2.0)
      sources.mgnify.max_file_mb:           int        (default: 2000)
    """

    BASE = "https://www.ebi.ac.uk/metagenomics/api/v1"

    def __init__(self, config, db):
        self.config = config
        self.db = db

        self.base_dir = Path(config['paths']['base_data']) / 'mgnify'
        self.base_dir.mkdir(parents=True, exist_ok=True)

        src_cfg = (config.get('sources', {}).get('mgnify', {}) or {})
        self.environments = src_cfg.get('environments', ['soil', 'marine', 'freshwater', 'plant', 'gut', 'sediment'])
        self.analyses_per_study = int(src_cfg.get('analyses_per_study', 2))
        self.delay = float(src_cfg.get('delay_seconds', 2.0))
        self.max_mb = int(src_cfg.get('max_file_mb', 2000))

        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "Mandrake-ETL-MGnify/1.0",
            "Accept": "application/json"
        })
        self.last_request = 0.0

        # analysis_id -> {'origin': ...}
        self.metadata = {}

    # ------------------- Public -------------------

    def download_batch(self, batch_size, seen_ids):
        downloaded = []
        for env in self.environments:
            if len(downloaded) >= batch_size:
                break

            studies = self._search_studies(env, limit=max(batch_size, 10))
            if not studies:
                continue

            for study in studies:
                if len(downloaded) >= batch_size:
                    break
                study_id = study.get('id')
                if not study_id:
                    continue

                analyses = self._list_analyses(study_id, page_size=self.analyses_per_study)
                if not analyses:
                    continue

                study_origin = self._infer_origin_from_study(study)

                for analysis in analyses:
                    if len(downloaded) >= batch_size:
                        break
                    analysis_id = analysis.get('id')
                    if not analysis_id or analysis_id in seen_ids or self.db.entry_exists(analysis_id):
                        continue

                    path, origin = self._download_fastq_fasta_only(analysis_id, study_origin)
                    if path:
                        downloaded.append(path)
                        seen_ids.add(analysis_id)
                        self.metadata[analysis_id] = {'origin': origin}

        logger.info(f"MGnify: Downloaded {len(downloaded)} metagenomes")
        return downloaded

    def get_metadata(self, accession_or_id: str):
        return self.metadata.get(accession_or_id, {})

    # ------------------- Internals -------------------

    def _rate(self):
        elapsed = time.time() - self.last_request
        if elapsed < self.delay:
            time.sleep(self.delay - elapsed)
        self.last_request = time.time()

    def _get_json(self, url, params=None, timeout=30):
        try:
            self._rate()
            r = self.session.get(url, params=params or {}, timeout=timeout)
            if r.ok:
                return r.json()
        except Exception:
            return None
        return None

    def _head(self, url, timeout=15):
        try:
            self._rate()
            return self.session.head(url, allow_redirects=True, timeout=timeout)
        except Exception:
            return None

    def _search_studies(self, environment: str, limit: int = 10):
        params = {"search": environment, "page_size": limit}
        data = self._get_json(f"{self.BASE}/studies", params=params)
        studies = (data or {}).get('data', [])
        logger.info(f"MGnify: found {len(studies)} studies for env='{environment}'")
        return studies

    def _list_analyses(self, study_id: str, page_size: int = 2):
        params = {"page_size": page_size}
        data = self._get_json(f"{self.BASE}/studies/{study_id}/analyses", params=params)
        return (data or {}).get('data', [])

    def _download_fastq_fasta_only(self, analysis_id: str, study_origin: str | None):
        """
        Accept ONLY items matching *FASTQ.fasta or *FASTQ.fasta.gz (case-insensitive).
        Skip everything else (.ffn, .faa, contigs).
        """
        dl = self._get_json(f"{self.BASE}/analyses/{analysis_id}/downloads")
        items = (dl or {}).get('data', [])
        if not items:
            return None, study_origin

        origin = self._infer_origin_from_analysis(analysis_id) or study_origin

        chosen = None
        for item in items:
            item_id = (item.get('id') or '')
            if not item_id:
                continue
            lid = item_id.lower()
            # exclude proteins/ORFs explicitly
            if lid.endswith('.faa') or '.faa' in lid or lid.endswith('.ffn') or '.ffn' in lid:
                continue
            # accept only processed reads exported as FASTA from FASTQ
            if ('fastq.fasta' in lid) or ('_fastq.fasta' in lid):
                chosen = item
                break

        if not chosen:
            # nothing acceptable for this analysis
            return None, origin

        url = (chosen.get('links') or {}).get('self')
        if not url:
            return None, origin

        # Build filename & dest folder
        filename = chosen.get('id')
        if not filename:
            return None, origin
        # Ensure extension consistency
        if not (filename.endswith('.fasta') or filename.endswith('.fasta.gz')):
            if filename.endswith('.gz'):
                filename = filename.replace('.gz', '.fasta.gz')
            else:
                filename += '.fasta'

        out_dir = self.base_dir / filename  # separate namespace per file root to avoid collisions in temp
        out_dir.mkdir(parents=True, exist_ok=True)
        dest = out_dir / filename

        path = self._download_stream(url, dest, max_size_mb=self.max_mb)
        if path:
            meta = {
                'analysis_id': analysis_id,
                'download_id': chosen.get('id', ''),
                'url': url,
                'origin': origin
            }
            with open(out_dir / f"{analysis_id}_download.json", 'w') as f:
                json.dump(meta, f, indent=2)
            return path, origin

        return None, origin

    def _download_stream(self, url: str, output_path: Path, max_size_mb: int):
        if output_path.exists() and output_path.stat().st_size > 0:
            logger.info(f"MGnify: exists {output_path.name}")
            return output_path

        # HEAD for size guard
        head = self._head(url)
        if head is not None:
            try:
                cl = int(head.headers.get('content-length', 0))
                if cl > 0 and (cl / (1024 * 1024)) > max_size_mb:
                    logger.warning(f"MGnify: skip {output_path.name} (>{max_size_mb} MB)")
                    return None
            except Exception:
                pass

        try:
            self._rate()
            with self.session.get(url, stream=True, timeout=300) as r:
                if not r.ok:
                    logger.warning(f"MGnify: download failed HTTP {r.status_code} for {url}")
                    return None

                tmp = output_path.with_suffix(output_path.suffix + ".tmp")
                with open(tmp, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            f.write(chunk)

                # Detect if payload is gz; keep gz if so
                try:
                    with gzip.open(tmp, 'rb') as gzf:
                        _ = gzf.read(2)
                    if not str(output_path).endswith('.gz'):
                        output_path = output_path.with_suffix(output_path.suffix + ".gz")
                except Exception:
                    pass

                tmp.rename(output_path)
                logger.info(f"MGnify: downloaded {output_path.name}")
                return output_path
        except requests.exceptions.Timeout:
            logger.warning(f"MGnify: timeout downloading {output_path.name}")
        except Exception as e:
            logger.warning(f"MGnify: download error {e} for {output_path.name}")
        return None

    # ------------------- Origin inference -------------------

    def _infer_origin_from_study(self, study_obj) -> str | None:
        attrs = (study_obj or {}).get('attributes', {}) or {}
        text = " ".join([
            str(attrs.get('biome') or ''),
            str(attrs.get('study-abstract') or ''),
            str(attrs.get('study-name') or ''),
        ]).lower()
        return self._normalize_origin(text)

    def _infer_origin_from_analysis(self, analysis_id: str) -> str | None:
        data = self._get_json(f"{self.BASE}/analyses/{analysis_id}")
        attrs = (data or {}).get('data', {}).get('attributes', {}) or {}
        text = " ".join([
            str(attrs.get('environment') or ''),
            str(attrs.get('environment_biome') or ''),
            str(attrs.get('biome') or ''),
            str(attrs.get('sample_desc') or ''),
            str(attrs.get('sample_name') or ''),
        ]).lower()
        return self._normalize_origin(text)

    @staticmethod
    def _normalize_origin(text: str) -> str | None:
        if not text:
            return None
        if any(k in text for k in ('soil',)):
            return 'soil'
        if any(k in text for k in ('marine', 'freshwater', 'water', 'ocean', 'sea', 'lake', 'river')):
            return 'water'
        if any(k in text for k in ('root', 'leaf', 'plant', 'rhizosphere')):
            return 'plant'
        if any(k in text for k in ('gut', 'fecal', 'feces', 'stool', 'intest')):
            return 'gut'
        if any(k in text for k in ('sediment', 'mud', 'silt')):
            return 'sediment'
        if any(k in text for k in ('skin', 'oral', 'mouth', 'saliva')):
            return 'host'
        return 'Unknown'
