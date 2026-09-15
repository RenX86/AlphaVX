"""Tests for the SQLite result cache module."""

import pytest
from pathlib import Path

from alphavx.cache import ResultCache


@pytest.fixture
def cache(tmp_path: Path) -> ResultCache:
    """Create a fresh cache in a temporary directory."""
    return ResultCache(tmp_path / "test_cache", config_hash="test_hash")


class TestResultCache:
    """Tests for ResultCache."""

    def test_put_and_get(self, cache: ResultCache) -> None:
        data = {"raw_score": [0.5, -0.3], "gene_name": ["TP53", "TP53"]}
        cache.put("chr17:7674220:G>A", data)

        result = cache.get("chr17:7674220:G>A")
        assert result is not None
        assert result["raw_score"] == [0.5, -0.3]
        assert result["gene_name"] == ["TP53", "TP53"]

    def test_get_missing_returns_none(self, cache: ResultCache) -> None:
        assert cache.get("chr1:100:A>G") is None

    def test_has(self, cache: ResultCache) -> None:
        assert not cache.has("chr1:100:A>G")
        cache.put("chr1:100:A>G", {"score": [1.0]})
        assert cache.has("chr1:100:A>G")

    def test_overwrite_existing(self, cache: ResultCache) -> None:
        cache.put("chr1:100:A>G", {"score": [1.0]})
        cache.put("chr1:100:A>G", {"score": [2.0]})

        result = cache.get("chr1:100:A>G")
        assert result["score"] == [2.0]

    def test_clear(self, cache: ResultCache) -> None:
        cache.put("chr1:100:A>G", {"score": [1.0]})
        cache.put("chr2:200:C>T", {"score": [2.0]})
        assert cache.stats()["count"] == 2

        cache.clear()
        assert cache.stats()["count"] == 0
        assert cache.get("chr1:100:A>G") is None

    def test_stats(self, cache: ResultCache) -> None:
        stats = cache.stats()
        assert stats["count"] == 0
        assert stats["cache_size_bytes"] >= 0

        cache.put("chr1:100:A>G", {"score": [1.0]})
        cache.put("chr2:200:C>T", {"score": [2.0]})

        stats = cache.stats()
        assert stats["count"] == 2
        assert stats["cache_size_bytes"] > 0

    def test_cache_dir_created(self, tmp_path: Path) -> None:
        cache_dir = tmp_path / "deep" / "nested" / "cache"
        assert not cache_dir.exists()

        cache = ResultCache(cache_dir)
        assert cache_dir.exists()
        assert cache.db_path.exists()


class TestCacheConfigInvalidation:
    """Tests for config-aware cache invalidation (bug fix #1)."""

    def test_different_config_hash_misses(self, tmp_path: Path) -> None:
        """Same variant with different config_hash should be a cache miss."""
        cache_dir = tmp_path / "cache"
        cache_a = ResultCache(cache_dir, config_hash="config_A")
        cache_a.put("chr1:100:A>G", {"score": [1.0]})

        cache_b = ResultCache(cache_dir, config_hash="config_B")
        assert not cache_b.has("chr1:100:A>G")
        assert cache_b.get("chr1:100:A>G") is None

    def test_same_config_hash_hits(self, tmp_path: Path) -> None:
        """Same variant with same config_hash should be a cache hit."""
        cache_dir = tmp_path / "cache"
        cache1 = ResultCache(cache_dir, config_hash="same_hash")
        cache1.put("chr1:100:A>G", {"score": [1.0]})

        cache2 = ResultCache(cache_dir, config_hash="same_hash")
        assert cache2.has("chr1:100:A>G")
        assert cache2.get("chr1:100:A>G") == {"score": [1.0]}

    def test_both_configs_coexist(self, tmp_path: Path) -> None:
        """Different configs should store independently for the same variant."""
        cache_dir = tmp_path / "cache"
        cache_a = ResultCache(cache_dir, config_hash="A")
        cache_b = ResultCache(cache_dir, config_hash="B")

        cache_a.put("chr1:100:A>G", {"score": [1.0]})
        cache_b.put("chr1:100:A>G", {"score": [2.0]})

        assert cache_a.get("chr1:100:A>G") == {"score": [1.0]}
        assert cache_b.get("chr1:100:A>G") == {"score": [2.0]}

    def test_none_config_hash_always_misses(self, tmp_path: Path) -> None:
        """Cache with config_hash=None should always miss (safe default)."""
        cache_dir = tmp_path / "cache"

        # Put something with a real hash
        cache_real = ResultCache(cache_dir, config_hash="real")
        cache_real.put("chr1:100:A>G", {"score": [1.0]})

        # None hash → miss
        cache_none = ResultCache(cache_dir, config_hash=None)
        assert not cache_none.has("chr1:100:A>G")
        assert cache_none.get("chr1:100:A>G") is None

    def test_legacy_db_migration(self, tmp_path: Path) -> None:
        """Opening a DB created without config_hash should auto-migrate."""
        import sqlite3

        cache_dir = tmp_path / "legacy_cache"
        cache_dir.mkdir()
        db_path = cache_dir / "alphavx_cache.db"

        # Create a legacy schema (no config_hash column)
        with sqlite3.connect(db_path) as conn:
            conn.execute("""
                CREATE TABLE variant_scores (
                    variant_key TEXT PRIMARY KEY,
                    scores_json TEXT NOT NULL,
                    scored_at TEXT NOT NULL
                )
            """)
            conn.execute(
                "INSERT INTO variant_scores VALUES (?, ?, ?)",
                ("chr1:100:A>G", '{"score": [1.0]}', "2024-01-01"),
            )
            conn.commit()

        # Opening with ResultCache should migrate without error
        cache = ResultCache(cache_dir, config_hash="new_hash")

        # Legacy rows have NULL config_hash → should be misses for any hash
        assert not cache.has("chr1:100:A>G")

        # But we can put and get new entries
        cache.put("chr1:100:A>G", {"score": [2.0]})
        assert cache.has("chr1:100:A>G")
        assert cache.get("chr1:100:A>G") == {"score": [2.0]}
