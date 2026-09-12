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
    variant_key TEXT PRIMARY KEY,
    scores_json TEXT NOT NULL,
    scored_at TEXT NOT NULL
);
"""


class ResultCache:
    """SQLite cache for variant scoring results.

    Args:
        cache_dir: Directory where the cache database will be stored.
    """

    def __init__(self, cache_dir: Path) -> None:
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.db_path = self.cache_dir / "alphavx_cache.db"
        self._init_db()

    def _init_db(self) -> None:
        """Create the cache table if it doesn't exist."""
        with sqlite3.connect(self.db_path) as conn:
            conn.execute(_SCHEMA)
            conn.commit()
        logger.debug("Cache initialized at %s", self.db_path)

    def has(self, variant_key: str) -> bool:
        """Check if a variant has cached results.

        Args:
            variant_key: Variant identifier (chr:pos:ref>alt).

        Returns:
            True if the variant has cached scores.
        """
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.execute(
                "SELECT 1 FROM variant_scores WHERE variant_key = ?",
                (variant_key,),
            )
            return cursor.fetchone() is not None

    def get(self, variant_key: str) -> dict | None:
        """Retrieve cached scores for a variant.

        Args:
            variant_key: Variant identifier (chr:pos:ref>alt).

        Returns:
            Deserialized scores dict, or None if not cached.
        """
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.execute(
                "SELECT scores_json FROM variant_scores WHERE variant_key = ?",
                (variant_key,),
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
        with sqlite3.connect(self.db_path) as conn:
            conn.execute(
                "INSERT OR REPLACE INTO variant_scores (variant_key, scores_json, scored_at) "
                "VALUES (?, ?, ?)",
                (variant_key, scores_json, now),
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
