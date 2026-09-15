"""Tests for the report generation module.

Covers:
- CSV report generation (full + significant)
- HTML report generation
- HTML escaping of user-controlled values (bug fix verification)
- Significance flagging
- Edge cases (empty DataFrame, missing columns)
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from alphavx.reporter import generate_csv_report, generate_html_report, generate_report


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_df(n_variants: int = 3, significant: bool = True) -> pd.DataFrame:
    """Build a realistic scoring DataFrame."""
    rows = []
    for i in range(n_variants):
        q_score = 0.998 if significant and i == 0 else 0.5
        rows.append({
            "variant_key": f"chr{i + 1}:{100 + i}:A>G",
            "gene_name": f"GENE{i}",
            "output_type": "RNA_SEQ",
            "biosample_name": f"tissue_{i}",
            "raw_score": 0.42 + i * 0.1,
            "quantile_score": q_score,
            "chrom": f"chr{i + 1}",
            "pos": 100 + i,
            "ref": "A",
            "alt": "G",
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# CSV report
# ---------------------------------------------------------------------------

class TestCsvReport:
    """Tests for generate_csv_report."""

    def test_creates_scores_csv(self, tmp_path: Path):
        df = _make_df()
        path = generate_csv_report(df, tmp_path)

        assert path.exists()
        assert path.name == "scores.csv"

        result = pd.read_csv(path)
        assert len(result) == 3
        assert "significant" in result.columns

    def test_creates_significant_csv(self, tmp_path: Path):
        df = _make_df(significant=True)
        generate_csv_report(df, tmp_path)

        sig_path = tmp_path / "significant.csv"
        assert sig_path.exists()

        sig_df = pd.read_csv(sig_path)
        assert len(sig_df) >= 1
        assert all(sig_df["significant"])

    def test_significance_threshold(self, tmp_path: Path):
        df = _make_df()
        generate_csv_report(df, tmp_path, quantile_threshold=0.99)

        sig_path = tmp_path / "significant.csv"
        sig_df = pd.read_csv(sig_path)
        # With threshold 0.99, the variant with q=0.998 should be significant
        assert len(sig_df) >= 1

    def test_no_quantile_column(self, tmp_path: Path):
        """Should handle missing quantile_score gracefully."""
        df = pd.DataFrame({"variant_key": ["v1"], "raw_score": [0.5]})
        path = generate_csv_report(df, tmp_path)

        result = pd.read_csv(path)
        assert not result["significant"].any()

    def test_creates_output_dir(self, tmp_path: Path):
        out = tmp_path / "deep" / "nested" / "results"
        df = _make_df()
        generate_csv_report(df, out)
        assert out.exists()


# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------

class TestHtmlReport:
    """Tests for generate_html_report."""

    def test_creates_report_html(self, tmp_path: Path):
        df = _make_df()
        path = generate_html_report(df, tmp_path)

        assert path.exists()
        assert path.name == "report.html"

        html = path.read_text(encoding="utf-8")
        assert "<!DOCTYPE html>" in html
        assert "AlphaVX" in html

    def test_contains_variant_count(self, tmp_path: Path):
        df = _make_df(n_variants=5)
        path = generate_html_report(df, tmp_path)
        html_content = path.read_text(encoding="utf-8")
        # The report renders total_variants inside a .value div
        assert ">5<" in html_content  # <div class="value">5</div>

    def test_html_escapes_values(self, tmp_path: Path):
        """Values containing HTML special chars must be escaped.

        This verifies the fix for the unescaped HTML injection bug.
        """
        df = pd.DataFrame({
            "variant_key": ['<script>alert("xss")</script>'],
            "gene_name": ["GENE&CO"],
            "output_type": ["RNA_SEQ"],
            "biosample_name": ['tissue<br>name'],
            "raw_score": [0.5],
            "quantile_score": [0.999],  # above threshold → in table
        })
        path = generate_html_report(df, tmp_path, quantile_threshold=0.995)
        html = path.read_text(encoding="utf-8")

        # Raw HTML should NOT appear unescaped
        assert "<script>" not in html
        assert "&lt;script&gt;" in html or "alert" not in html
        # Ampersand in gene name should be escaped
        assert "GENE&amp;CO" in html or "GENE&CO" not in html.replace("&amp;", "")

    def test_html_escapes_title(self, tmp_path: Path):
        """Title parameter should be escaped in the HTML output."""
        df = _make_df()
        path = generate_html_report(df, tmp_path, title='<img src=x onerror="alert(1)">')
        html = path.read_text(encoding="utf-8")

        assert 'onerror="alert(1)"' not in html

    def test_empty_dataframe(self, tmp_path: Path):
        df = pd.DataFrame()
        path = generate_html_report(df, tmp_path)
        assert path.exists()
        html = path.read_text(encoding="utf-8")
        assert "No significant results" in html


# ---------------------------------------------------------------------------
# generate_report (from scores.csv)
# ---------------------------------------------------------------------------

class TestGenerateReport:
    """Tests for generate_report (reads existing scores.csv)."""

    def test_generates_from_existing_scores(self, tmp_path: Path):
        # First create a scores.csv
        df = _make_df()
        df.to_csv(tmp_path / "scores.csv", index=False)

        generate_report(tmp_path)

        assert (tmp_path / "report.html").exists()

    def test_missing_scores_raises(self, tmp_path: Path):
        with pytest.raises(FileNotFoundError, match="No scores.csv"):
            generate_report(tmp_path)
