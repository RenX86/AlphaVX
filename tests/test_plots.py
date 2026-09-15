"""Tests for the visualization / plotting module.

Covers:
- Summary heatmap generation
- Per-variant detail plots
- Empty data / no significant variants placeholders
- Missing columns edge cases
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from alphavx.plots import plot_summary_heatmap, plot_variant_detail


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_df(n_variants: int = 3, significant: bool = True) -> pd.DataFrame:
    """Build a DataFrame with multiple variants and modalities."""
    rows = []
    modalities = ["RNA_SEQ", "DNASE", "ATAC"]
    for i in range(n_variants):
        for mod in modalities:
            q = 0.998 if significant and i == 0 else 0.3
            rows.append({
                "variant_key": f"chr{i + 1}:{100 + i}:A>G",
                "gene_name": f"GENE{i}",
                "output_type": mod,
                "biosample_name": f"tissue_{i}",
                "raw_score": 0.42 + i * 0.1,
                "quantile_score": q,
                "chrom": f"chr{i + 1}",
                "pos": 100 + i,
                "ref": "A",
                "alt": "G",
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Summary heatmap
# ---------------------------------------------------------------------------

class TestSummaryHeatmap:
    """Tests for plot_summary_heatmap."""

    def test_creates_png(self, tmp_path: Path):
        df = _make_df()
        out = tmp_path / "plots" / "heatmap.png"
        result = plot_summary_heatmap(df, out)

        assert result == out
        assert out.exists()
        assert out.stat().st_size > 0

    def test_no_significant_variants(self, tmp_path: Path):
        """Should create a placeholder when no variants pass threshold."""
        df = _make_df(significant=False)
        out = tmp_path / "heatmap.png"
        result = plot_summary_heatmap(df, out, quantile_threshold=0.999)

        assert out.exists()

    def test_empty_dataframe(self, tmp_path: Path):
        out = tmp_path / "heatmap.png"
        result = plot_summary_heatmap(pd.DataFrame(), out)
        assert out.exists()

    def test_missing_columns(self, tmp_path: Path):
        """Should handle DataFrame without expected columns."""
        df = pd.DataFrame({"x": [1, 2, 3]})
        out = tmp_path / "heatmap.png"
        result = plot_summary_heatmap(df, out)
        assert out.exists()

    def test_creates_parent_dirs(self, tmp_path: Path):
        out = tmp_path / "deep" / "nested" / "heatmap.png"
        plot_summary_heatmap(_make_df(), out)
        assert out.exists()


# ---------------------------------------------------------------------------
# Per-variant detail plot
# ---------------------------------------------------------------------------

class TestVariantDetail:
    """Tests for plot_variant_detail."""

    def test_creates_png(self, tmp_path: Path):
        df = _make_df()
        vk = "chr1:100:A>G"
        out = tmp_path / "variant.png"
        result = plot_variant_detail(df, vk, out)

        assert result == out
        assert out.exists()
        assert out.stat().st_size > 0

    def test_nonexistent_variant(self, tmp_path: Path):
        """Should produce a placeholder for unknown variant key."""
        df = _make_df()
        out = tmp_path / "missing.png"
        result = plot_variant_detail(df, "chrX:999:C>T", out)
        assert out.exists()

    def test_significance_coloring(self, tmp_path: Path):
        """Should not crash when mixing significant/non-significant bars."""
        df = _make_df(significant=True)
        out = tmp_path / "mixed.png"
        plot_variant_detail(df, "chr1:100:A>G", out, quantile_threshold=0.995)
        assert out.exists()

    def test_creates_parent_dirs(self, tmp_path: Path):
        df = _make_df()
        out = tmp_path / "a" / "b" / "plot.png"
        plot_variant_detail(df, "chr1:100:A>G", out)
        assert out.exists()
