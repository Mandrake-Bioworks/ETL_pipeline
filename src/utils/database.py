#!/usr/bin/env python3
"""Database Manager with filtering statistics"""
import logging
from datetime import datetime
from decimal import Decimal

from psycopg2.pool import SimpleConnectionPool
from psycopg2 import errors

logger = logging.getLogger(__name__)


class Database:
    def __init__(self, config: dict):
        dsn = config["aws"]["rds"]["connection_string"]
        max_conns = int(config.get("processing", {}).get("db_max_connections", 12))
        self.pool = SimpleConnectionPool(minconn=1, maxconn=max_conns, dsn=dsn)
        self._ensure_schema()
        logger.info(f"Database connected (pool size 1..{max_conns})")

    def _getconn(self):
        return self.pool.getconn()

    def _putconn(self, conn):
        self.pool.putconn(conn)

    def _ensure_schema(self):
        conn = self._getconn()
        try:
            prev_ac = conn.autocommit
            conn.autocommit = True
            with conn.cursor() as cur:
                # Main entries table
                cur.execute("""
                    CREATE TABLE IF NOT EXISTS entries (
                        id SERIAL PRIMARY KEY,
                        source VARCHAR(50),
                        accession VARCHAR(100),
                        s3_genome_path TEXT,
                        s3_protein_path TEXT,
                        sequence_hash VARCHAR(64),
                        total_bp BIGINT,
                        species TEXT,
                        kingdom TEXT,
                        origin TEXT,
                        status VARCHAR(50),
                        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                    )
                """)
                
                # Filtering statistics table
                cur.execute("""
                    CREATE TABLE IF NOT EXISTS filtering_stats (
                        id SERIAL PRIMARY KEY,
                        source VARCHAR(50),
                        accession VARCHAR(100),
                        total_contigs INT,
                        contigs_kept INT,
                        contigs_removed INT,
                        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                    )
                """)
                
                # Add columns if they don't exist
                cur.execute("ALTER TABLE entries ADD COLUMN IF NOT EXISTS source VARCHAR(50)")
                cur.execute("ALTER TABLE entries ADD COLUMN IF NOT EXISTS accession VARCHAR(100)")
                cur.execute("ALTER TABLE entries ADD COLUMN IF NOT EXISTS s3_genome_path TEXT")
                cur.execute("ALTER TABLE entries ADD COLUMN IF NOT EXISTS s3_protein_path TEXT")
                cur.execute("ALTER TABLE entries ADD COLUMN IF NOT EXISTS sequence_hash VARCHAR(64)")
                cur.execute("ALTER TABLE entries ADD COLUMN IF NOT EXISTS total_bp BIGINT")
                cur.execute("ALTER TABLE entries ADD COLUMN IF NOT EXISTS species TEXT")
                cur.execute("ALTER TABLE entries ADD COLUMN IF NOT EXISTS kingdom TEXT")
                cur.execute("ALTER TABLE entries ADD COLUMN IF NOT EXISTS origin TEXT")
                cur.execute("ALTER TABLE entries ADD COLUMN IF NOT EXISTS status VARCHAR(50)")
                cur.execute("ALTER TABLE entries ADD COLUMN IF NOT EXISTS created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP")

                # Indexes
                cur.execute("CREATE INDEX IF NOT EXISTS idx_entries_accession ON entries(accession)")
                cur.execute("CREATE INDEX IF NOT EXISTS idx_entries_species ON entries(species)")
                cur.execute("CREATE INDEX IF NOT EXISTS idx_entries_kingdom ON entries(kingdom)")
                cur.execute("CREATE INDEX IF NOT EXISTS idx_entries_origin ON entries(origin)")
                cur.execute("CREATE INDEX IF NOT EXISTS idx_entries_source ON entries(source)")

                # Unique constraints
                try:
                    cur.execute("CREATE UNIQUE INDEX IF NOT EXISTS uq_entries_sequence_hash ON entries(sequence_hash)")
                except Exception:
                    pass
                
                try:
                    cur.execute("CREATE UNIQUE INDEX IF NOT EXISTS uq_entries_source_accession ON entries(source, accession)")
                except errors.UniqueViolation:
                    logger.warning("Duplicate (source, accession) rows found; auto-deduplicating...")
                    cur.execute("""
                        DELETE FROM entries e
                        USING (
                            SELECT source, accession, MIN(id) AS keep_id
                            FROM entries
                            WHERE source IS NOT NULL AND accession IS NOT NULL
                            GROUP BY source, accession
                            HAVING COUNT(*) > 1
                        ) d
                        WHERE e.source = d.source
                          AND e.accession = d.accession
                          AND e.id <> d.keep_id
                    """)
                    cur.execute("CREATE UNIQUE INDEX IF NOT EXISTS uq_entries_source_accession ON entries(source, accession)")
                    logger.info("Dedup complete and unique (source, accession) index created.")

                # Source state table
                cur.execute("""
                    CREATE TABLE IF NOT EXISTS source_state (
                        source VARCHAR(50) PRIMARY KEY,
                        total_entries INT DEFAULT 0,
                        total_bp BIGINT DEFAULT 0,
                        last_update TIMESTAMP
                    )
                """)
            conn.autocommit = prev_ac
            logger.info("Schema ensured with filtering_stats table")
        finally:
            self._putconn(conn)

    def insert_entry(self, *, source, accession, s3_genome_path, s3_protein_path,
                     sequence_hash, total_bp, species=None, kingdom=None, origin=None) -> tuple[bool, str | None]:
        conn = self._getconn()
        try:
            with conn.cursor() as cur:
                cur.execute(
                    """
                    INSERT INTO entries (
                        source, accession, s3_genome_path, s3_protein_path,
                        sequence_hash, total_bp, species, kingdom, origin, status
                    )
                    VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, 'uploaded')
                    ON CONFLICT DO NOTHING
                    RETURNING id
                    """,
                    (source, accession, s3_genome_path, s3_protein_path,
                     sequence_hash, total_bp, species, kingdom, origin)
                )
                row = cur.fetchone()
                conn.commit()
                if row:
                    return True, None

                reason = None
                if sequence_hash is not None:
                    cur.execute("SELECT 1 FROM entries WHERE sequence_hash=%s LIMIT 1", (sequence_hash,))
                    if cur.fetchone():
                        reason = "hash_conflict"
                if reason is None:
                    cur.execute("SELECT 1 FROM entries WHERE source=%s AND accession=%s LIMIT 1", (source, accession))
                    if cur.fetchone():
                        reason = "accession_conflict"
                return False, (reason or "conflict")
        finally:
            self._putconn(conn)

    def insert_filtering_stats(self, source: str, accession: str, total_contigs: int, 
                               contigs_kept: int, contigs_removed: int):
        """Record filtering statistics"""
        conn = self._getconn()
        try:
            with conn.cursor() as cur:
                cur.execute(
                    """
                    INSERT INTO filtering_stats (source, accession, total_contigs, contigs_kept, contigs_removed)
                    VALUES (%s, %s, %s, %s, %s)
                    """,
                    (source, accession, total_contigs, contigs_kept, contigs_removed)
                )
                conn.commit()
        finally:
            self._putconn(conn)

    def entry_exists(self, accession: str) -> bool:
        conn = self._getconn()
        try:
            with conn.cursor() as cur:
                cur.execute("SELECT 1 FROM entries WHERE accession=%s LIMIT 1", (accession,))
                return cur.fetchone() is not None
        finally:
            self._putconn(conn)

    def get_existing_hashes(self):
        conn = self._getconn()
        try:
            with conn.cursor() as cur:
                cur.execute("SELECT sequence_hash FROM entries WHERE sequence_hash IS NOT NULL")
                return {row[0] for row in cur.fetchall()}
        finally:
            self._putconn(conn)

    def get_stats(self):
        conn = self._getconn()
        try:
            with conn.cursor() as cur:
                cur.execute("""
                    SELECT source, COUNT(*) AS entries,
                           COALESCE(SUM(total_bp), 0) AS bp,
                           COUNT(DISTINCT species) AS species_count
                    FROM entries
                    GROUP BY source
                    ORDER BY source
                """)
                rows = cur.fetchall()
            out = []
            for source, entries, bp, species in rows:
                bp_int = int(bp) if isinstance(bp, Decimal) else (bp or 0)
                out.append((source, entries, bp_int, species))
            return out
        finally:
            self._putconn(conn)

    def get_counts_by_kingdom(self):
        conn = self._getconn()
        try:
            with conn.cursor() as cur:
                cur.execute("""
                    SELECT COALESCE(NULLIF(TRIM(kingdom), ''), 'Unknown') AS k,
                           COUNT(*) AS entries,
                           COALESCE(SUM(total_bp), 0) AS bp,
                           COUNT(DISTINCT species) AS species_count
                    FROM entries
                    WHERE source IN ('ncbi','ena')
                    GROUP BY COALESCE(NULLIF(TRIM(kingdom), ''), 'Unknown')
                    ORDER BY entries DESC
                """)
                rows = cur.fetchall()
            return [(r[0], r[1], int(r[2]) if isinstance(r[2], Decimal) else (r[2] or 0), r[3]) for r in rows]
        finally:
            self._putconn(conn)

    def get_counts_by_origin(self):
        conn = self._getconn()
        try:
            with conn.cursor() as cur:
                cur.execute("""
                    SELECT COALESCE(NULLIF(TRIM(origin), ''), 'Unknown') AS o,
                           COUNT(*) AS entries,
                           COALESCE(SUM(total_bp), 0) AS bp,
                           COUNT(DISTINCT species) AS species_count
                    FROM entries
                    WHERE source = 'mgnify'
                    GROUP BY COALESCE(NULLIF(TRIM(origin), ''), 'Unknown')
                    ORDER BY entries DESC
                """)
                rows = cur.fetchall()
            return [(r[0], r[1], int(r[2]) if isinstance(r[2], Decimal) else (r[2] or 0), r[3]) for r in rows]
        finally:
            self._putconn(conn)

    def get_dedup_stats(self):
        """Get simplified deduplication statistics"""
        conn = self._getconn()
        try:
            with conn.cursor() as cur:
                # Overall stats
                cur.execute("""
                    SELECT 
                        COUNT(*) AS total_entries,
                        COUNT(DISTINCT COALESCE(sequence_hash, accession)) AS unique_entries,
                        COUNT(*) - COUNT(DISTINCT COALESCE(sequence_hash, accession)) AS duplicate_entries
                    FROM entries
                """)
                total_entries, unique_entries, duplicate_entries = cur.fetchone()
            
            return {
                'total_entries': total_entries or 0,
                'unique_entries': unique_entries or 0,
                'duplicate_entries': duplicate_entries or 0
            }
        finally:
            self._putconn(conn)

    def get_filtering_stats(self):
        """Get filtering statistics"""
        conn = self._getconn()
        try:
            with conn.cursor() as cur:
                cur.execute("""
                    SELECT 
                        COALESCE(SUM(total_contigs), 0) AS total_contigs,
                        COALESCE(SUM(contigs_kept), 0) AS contigs_kept,
                        COALESCE(SUM(contigs_removed), 0) AS contigs_removed
                    FROM filtering_stats
                """)
                row = cur.fetchone()
                if row:
                    return {
                        'total_contigs': int(row[0]),
                        'contigs_kept': int(row[1]),
                        'contigs_removed': int(row[2])
                    }
            return {'total_contigs': 0, 'contigs_kept': 0, 'contigs_removed': 0}
        finally:
            self._putconn(conn)

    def close(self):
        try:
            self.pool.closeall()
        except Exception:
            pass