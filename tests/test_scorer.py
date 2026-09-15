"""Tests for the variant scoring engine.

Covers:
- Single-variant scoring with mocked AlphaGenome client
- Batch scoring loop
- Cache hit / miss / store behaviour
- Retry logic with exponential backoff
- DataFrame assembly and column contract
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch, PropertyMock
import time

import pandas as pd
import pytest

from alphavx.config import Config
from alphavx.scorer import VariantScorer
from alphavx.vcf_parser import VariantRecord


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_record(chrom="chr17", pos=7674220, ref="G", alt="A") -> VariantRecord:
    return VariantRecord(chrom=chrom, pos=pos, ref=ref, alt=alt)


def _make_config(**overrides) -> Config:
    defaults = dict(
        api_key="test-key",
        max_retries=2,
        retry_delay=0.01,  # keep tests fast
        sequence_length=1024,
        quantile_threshold=0.995,
        modalities=["RNA_SEQ", "DNASE"],
    )
    defaults.update(overrides)
    return Config(**defaults)


def _mock_score_adata():
    """Return a minimal mock that variant_scorers.tidy_scores can process."""
    return MagicMock()


def _fake_tidy_df(variant_key: str = "chr17:7674220:G>A") -> pd.DataFrame:
    """Build a realistic tidy DataFrame as score_variant would produce."""
    return pd.DataFrame({
        "biosample_name": ["brain", "liver"],
        "gene_name": ["TP53", "TP53"],
        "output_type": ["RNA_SEQ", "DNASE"],
        "raw_score": [0.42, -0.18],
        "quantile_score": [0.997, 0.500],
    })


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def config() -> Config:
    return _make_config()


@pytest.fixture
def scorer(config: Config) -> VariantScorer:
    return VariantScorer(config)


@pytest.fixture
def mock_cache():
    """A MagicMock that behaves like ResultCache."""
    cache = MagicMock()
    cache.has.return_value = False
    cache.get.return_value = None
    return cache


# ---------------------------------------------------------------------------
# VariantScorer — score_variant
# ---------------------------------------------------------------------------

class TestScoreVariant:
    """Tests for VariantScorer.score_variant."""

    @patch("alphavx.scorer.VariantScorer._init_client")
    @patch("alphavx.scorer.VariantScorer._get_scorers", return_value=["mock_scorer"])
    def test_returns_dataframe_with_variant_columns(self, _scorers, _init, scorer):
        """score_variant should return a DF with chrom/pos/ref/alt/variant_key."""
        import sys
        record = _make_record()

        # Mock the client.score_variant call
        mock_adata = _mock_score_adata()
        scorer._client = MagicMock()
        scorer._client.score_variant.return_value = [mock_adata]

        # Set up the mock module attributes that score_variant imports lazily
        genome_mock = sys.modules["alphagenome.data.genome"]
        variant_scorers_mock = sys.modules["alphagenome.models.variant_scorers"]
        variant_scorers_mock.tidy_scores = MagicMock(return_value=_fake_tidy_df())

        df = scorer.score_variant(record)

        assert not df.empty
        for col in ("chrom", "pos", "ref", "alt", "variant_key"):
            assert col in df.columns, f"Missing column: {col}"
        assert (df["chrom"] == "chr17").all()
        assert (df["pos"] == 7674220).all()
        assert (df["variant_key"] == record.key).all()

    @patch("alphavx.scorer.VariantScorer._init_client")
    @patch("alphavx.scorer.VariantScorer._get_scorers", return_value=["mock_scorer"])
    def test_returns_empty_df_when_no_scores(self, _scorers, _init, scorer):
        """score_variant returns empty DF when API returns nothing."""
        record = _make_record()
        scorer._client = MagicMock()
        scorer._client.score_variant.return_value = []

        df = scorer.score_variant(record)

        assert df.empty


# ---------------------------------------------------------------------------
# VariantScorer — score_batch
# ---------------------------------------------------------------------------

class TestScoreBatch:
    """Tests for VariantScorer.score_batch (loop, cache, progress)."""

    def test_batch_returns_combined_df(self, scorer, mock_cache):
        """score_batch should combine results from all variants."""
        records = [_make_record(pos=100 + i) for i in range(3)]
        fake_dfs = [
            pd.DataFrame({"raw_score": [float(i)], "variant_key": [r.key]})
            for i, r in enumerate(records)
        ]

        with patch.object(scorer, "_score_with_retry", side_effect=fake_dfs):
            df = scorer.score_batch(records, cache=mock_cache)

        assert len(df) == 3
        assert set(df["variant_key"]) == {r.key for r in records}

    def test_batch_uses_cache_hits(self, scorer, mock_cache):
        """Cached variants should be loaded from cache, not scored."""
        record = _make_record()
        cached_data = {"raw_score": [0.99], "variant_key": [record.key]}

        mock_cache.has.return_value = True
        mock_cache.get.return_value = cached_data

        with patch.object(scorer, "_score_with_retry") as mock_score:
            df = scorer.score_batch([record], cache=mock_cache)

        # Should NOT have called the scorer
        mock_score.assert_not_called()
        assert len(df) == 1
        assert df["raw_score"].iloc[0] == 0.99

    def test_batch_stores_results_in_cache(self, scorer, mock_cache):
        """Successfully scored variants should be cached."""
        record = _make_record()
        fake_df = pd.DataFrame({"raw_score": [0.5], "variant_key": [record.key]})

        with patch.object(scorer, "_score_with_retry", return_value=fake_df):
            scorer.score_batch([record], cache=mock_cache)

        mock_cache.put.assert_called_once()
        call_key = mock_cache.put.call_args[0][0]
        assert call_key == record.key

    def test_batch_continues_on_error(self, scorer, mock_cache):
        """Errors for individual variants shouldn't stop the batch."""
        records = [_make_record(pos=100), _make_record(pos=200)]
        fake_results = [None, pd.DataFrame({"raw_score": [0.5], "variant_key": [records[1].key]})]

        with patch.object(scorer, "_score_with_retry", side_effect=fake_results):
            df = scorer.score_batch(records, cache=mock_cache)

        assert len(df) == 1
        assert df["variant_key"].iloc[0] == records[1].key

    def test_batch_returns_empty_df_when_all_fail(self, scorer, mock_cache):
        """If all variants fail, return empty DataFrame."""
        records = [_make_record(pos=100), _make_record(pos=200)]

        with patch.object(scorer, "_score_with_retry", return_value=None):
            df = scorer.score_batch(records, cache=mock_cache)

        assert df.empty

    def test_batch_progress_callback(self, scorer, mock_cache):
        """Progress callback should be invoked for each variant."""
        records = [_make_record(pos=i) for i in range(5)]
        calls = []

        def track(i, total, rec):
            calls.append((i, total, rec.key))

        with patch.object(scorer, "_score_with_retry", return_value=None):
            scorer.score_batch(records, cache=mock_cache, progress_callback=track)

        assert len(calls) == 5
        assert all(total == 5 for _, total, _ in calls)
        assert [i for i, _, _ in calls] == [0, 1, 2, 3, 4]

    def test_batch_without_cache(self, scorer):
        """score_batch should work fine when cache is None."""
        record = _make_record()
        fake_df = pd.DataFrame({"raw_score": [0.5], "variant_key": [record.key]})

        with patch.object(scorer, "_score_with_retry", return_value=fake_df):
            df = scorer.score_batch([record], cache=None)

        assert len(df) == 1


# ---------------------------------------------------------------------------
# VariantScorer — retry logic
# ---------------------------------------------------------------------------

class TestRetryLogic:
    """Tests for _score_with_retry."""

    def test_returns_on_first_success(self, scorer):
        """No retry needed when first attempt succeeds."""
        record = _make_record()
        expected = pd.DataFrame({"raw_score": [1.0]})

        with patch.object(scorer, "score_variant", return_value=expected) as mock_sv:
            result = scorer._score_with_retry(record)

        mock_sv.assert_called_once()
        assert result is not None
        assert result["raw_score"].iloc[0] == 1.0

    def test_retries_on_failure_then_succeeds(self, scorer):
        """Should retry and return result when a later attempt succeeds."""
        record = _make_record()
        expected = pd.DataFrame({"raw_score": [1.0]})

        with patch.object(
            scorer, "score_variant",
            side_effect=[RuntimeError("fail"), RuntimeError("fail"), expected],
        ) as mock_sv:
            result = scorer._score_with_retry(record)

        assert mock_sv.call_count == 3  # 1 initial + 2 retries (max_retries=2)
        assert result is not None

    def test_returns_none_after_all_retries_exhausted(self, scorer):
        """Should return None when all attempts fail."""
        record = _make_record()

        with patch.object(
            scorer, "score_variant",
            side_effect=RuntimeError("persistent failure"),
        ) as mock_sv:
            result = scorer._score_with_retry(record)

        assert mock_sv.call_count == 3  # 1 + max_retries(2)
        assert result is None

    def test_retry_respects_max_retries_config(self):
        """Number of attempts = max_retries + 1."""
        config = _make_config(max_retries=5, retry_delay=0.001)
        scorer = VariantScorer(config)
        record = _make_record()

        with patch.object(
            scorer, "score_variant",
            side_effect=RuntimeError("fail"),
        ) as mock_sv:
            scorer._score_with_retry(record)

        assert mock_sv.call_count == 6  # 1 + 5


# ---------------------------------------------------------------------------
# VariantScorer — _get_scorers
# ---------------------------------------------------------------------------

class TestGetScorers:
    """Tests for _get_scorers modality lookup."""

    def test_known_modalities_returned(self):
        import sys
        config = _make_config(modalities=["RNA_SEQ", "DNASE"])
        scorer = VariantScorer(config)

        mock_scorers = {"RNA_SEQ": "rna_scorer", "DNASE": "dnase_scorer", "ATAC": "atac_scorer"}
        vs_mod = sys.modules["alphagenome.models.variant_scorers"]
        original = vs_mod.RECOMMENDED_VARIANT_SCORERS
        try:
            vs_mod.RECOMMENDED_VARIANT_SCORERS = mock_scorers
            result = scorer._get_scorers()
        finally:
            vs_mod.RECOMMENDED_VARIANT_SCORERS = original

        assert result == ["rna_scorer", "dnase_scorer"]

    def test_unknown_modality_skipped(self):
        import sys
        config = _make_config(modalities=["RNA_SEQ", "NONEXISTENT"])
        scorer = VariantScorer(config)

        mock_scorers = {"RNA_SEQ": "rna_scorer"}
        vs_mod = sys.modules["alphagenome.models.variant_scorers"]
        original = vs_mod.RECOMMENDED_VARIANT_SCORERS
        try:
            vs_mod.RECOMMENDED_VARIANT_SCORERS = mock_scorers
            result = scorer._get_scorers()
        finally:
            vs_mod.RECOMMENDED_VARIANT_SCORERS = original

        assert result == ["rna_scorer"]

