"""SQLite-based result caching for AlphaVX.

Caches variant scoring results so interrupted batch runs can resume
without re-querying the AlphaGenome API for already-scored variants.
"""

from __future__ import annotations

import json
import logging
import os
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

logger = logging.getLogger(__name__)

_SCHEMA = """
CREATE TABLE IF NOT EXISTS variant_scores (
    variant_key TEXT NOT NULL,
    scores_json TEXT NOT NULL,
    scored_at TEXT NOT NULL,
    config_hash TEXT
);
"""

# Index for fast lookups on the composite key used by has() / get().
_INDEX = """
CREATE INDEX IF NOT EXISTS idx_variant_config
ON variant_scores (variant_key, config_hash);
"""


class ResultCache:
    """SQLite cache for variant scoring results.

    Results are keyed by both the variant identifier **and** an optional
    config hash so that different scoring configurations (modalities,
    sequence_length, etc.) are cached independently.

    Args:
        cache_dir: Directory where the cache database will be stored.
        config_hash: Optional hash string identifying the scoring
            configuration.  When *None*, every ``has``/``get`` call will
            behave as a cache miss (safe default for backward compat).
    """

    def __init__(self, cache_dir: Path, config_hash: str | None = None) -> None:
        import threading
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.db_path = self.cache_dir / "alphavx_cache.db"
        self.config_hash = config_hash
        self._lock = threading.Lock()
        self._init_db()

    def _init_db(self) -> None:
        """Create the cache table if it doesn't exist and migrate old schemas."""
        with sqlite3.connect(self.db_path) as conn:
            conn.execute(_SCHEMA)
            # Migrate legacy databases that lack the config_hash column.
            self._migrate(conn)
            conn.execute(_INDEX)
            conn.commit()
        logger.debug("Cache initialized at %s", self.db_path)

    @staticmethod
    def _migrate(conn: sqlite3.Connection) -> None:
        """Add the config_hash column if it is missing (legacy DB migration)."""
        cursor = conn.execute("PRAGMA table_info(variant_scores)")
        columns = {row[1] for row in cursor.fetchall()}
        if "config_hash" not in columns:
            # Recreate table to remove the PRIMARY KEY constraint on variant_key
            conn.execute("ALTER TABLE variant_scores RENAME TO variant_scores_old")
            conn.execute(_SCHEMA)
            conn.execute(
                "INSERT INTO variant_scores (variant_key, scores_json, scored_at, config_hash) "
                "SELECT variant_key, scores_json, scored_at, NULL FROM variant_scores_old"
            )
            conn.execute("DROP TABLE variant_scores_old")
            logger.info("Migrated cache DB: added config_hash column and removed old primary key")

    def has(self, variant_key: str) -> bool:
        """Check if a variant has cached results for the current config.

        Args:
            variant_key: Variant identifier (chr:pos:ref>alt).

        Returns:
            True if the variant has cached scores **and** the cache was
            created with a non-None ``config_hash``.
        """
        if self.config_hash is None:
            return False
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.execute(
                "SELECT 1 FROM variant_scores "
                "WHERE variant_key = ? AND config_hash = ?",
                (variant_key, self.config_hash),
            )
            return cursor.fetchone() is not None

    def get(self, variant_key: str) -> dict | None:
        """Retrieve cached scores for a variant under the current config.

        Args:
            variant_key: Variant identifier (chr:pos:ref>alt).

        Returns:
            Deserialized scores dict, or None if not cached.
        """
        if self.config_hash is None:
            return None
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.execute(
                "SELECT scores_json FROM variant_scores "
                "WHERE variant_key = ? AND config_hash = ?",
                (variant_key, self.config_hash),
            )
            row = cursor.fetchone()
            if row is None:
                return None
            return json.loads(row[0])

    def put(self, variant_key: str, scores: dict) -> None:
        """Store scores for a variant in the cache.

        Args:
            variant_key: Variant identifier (chr:pos:ref>alt).
            scores: Scoring results to cache (must be JSON-serializable).
        """
        now = datetime.now(timezone.utc).isoformat()
        scores_json = json.dumps(scores, default=str)
        with self._lock:
            with sqlite3.connect(self.db_path) as conn:
                # Remove any previous entry for the same (variant, config) pair.
                conn.execute(
                    "DELETE FROM variant_scores "
                    "WHERE variant_key = ? AND config_hash IS ?",
                    (variant_key, self.config_hash),
                )
                conn.execute(
                    "INSERT INTO variant_scores "
                    "(variant_key, scores_json, scored_at, config_hash) "
                    "VALUES (?, ?, ?, ?)",
                    (variant_key, scores_json, now, self.config_hash),
                )
            conn.commit()
        logger.debug("Cached scores for %s", variant_key)

    def clear(self) -> None:
        """Delete all cached results."""
        with sqlite3.connect(self.db_path) as conn:
            conn.execute("DELETE FROM variant_scores")
            conn.commit()
        logger.info("Cache cleared")

    def stats(self) -> dict:
        """Get cache statistics.

        Returns:
            Dict with 'count' (number of cached variants) and
            'cache_size_bytes' (database file size).
        """
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.execute("SELECT COUNT(*) FROM variant_scores")
            count = cursor.fetchone()[0]

        size = os.path.getsize(self.db_path) if self.db_path.exists() else 0
        return {"count": count, "cache_size_bytes": size}
