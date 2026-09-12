"""Tests for the SQLite result cache module."""

import pytest
from pathlib import Path

from alphavx.cache import ResultCache


@pytest.fixture
def cache(tmp_path: Path) -> ResultCache:
    """Create a fresh cache in a temporary directory."""
    return ResultCache(tmp_path / "test_cache")


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
