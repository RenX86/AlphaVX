"""Report generation for AlphaVX — CSV and self-contained HTML output."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def generate_csv_report(
    df: pd.DataFrame,
    output_dir: Path,
    quantile_threshold: float = 0.995,
) -> Path:
    """Save full and significant-only CSV reports.

    Args:
        df: Complete scoring results DataFrame.
        output_dir: Directory to write output files.
        quantile_threshold: Cutoff for significance flagging.

    Returns:
        Path to the full scores.csv file.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Add significance flag
    if "quantile_score" in df.columns:
        df["significant"] = df["quantile_score"].abs() > quantile_threshold
    else:
        df["significant"] = False

    # Full results
    scores_path = output_dir / "scores.csv"
    df.to_csv(scores_path, index=False)
    logger.info("Full results: %s (%d rows)", scores_path, len(df))

    # Significant only
    sig_df = df[df["significant"]]
    sig_path = output_dir / "significant.csv"
    sig_df.to_csv(sig_path, index=False)
    logger.info("Significant results: %s (%d rows)", sig_path, len(sig_df))

    return scores_path


def generate_html_report(
    df: pd.DataFrame,
    output_dir: Path,
    title: str = "AlphaVX Variant Report",
    quantile_threshold: float = 0.995,
) -> Path:
    """Generate a self-contained HTML report from scoring results.

    Args:
        df: Complete scoring results DataFrame.
        output_dir: Directory to write the report.
        title: Report title.
        quantile_threshold: Cutoff for significance highlighting.

    Returns:
        Path to the generated report.html file.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Compute summary stats
    if "quantile_score" in df.columns:
        sig_mask = df["quantile_score"].abs() > quantile_threshold
    else:
        sig_mask = pd.Series([False] * len(df))

    total_variants = df["variant_key"].nunique() if "variant_key" in df.columns else 0
    sig_count = df.loc[sig_mask, "variant_key"].nunique() if "variant_key" in df.columns else 0
    modalities = sorted(df["output_type"].unique().tolist()) if "output_type" in df.columns else []

    # Build significant hits table rows
    sig_df = df[sig_mask].copy() if sig_mask.any() else pd.DataFrame()
    if not sig_df.empty and "quantile_score" in sig_df.columns:
        sig_df = sig_df.sort_values("quantile_score", key=abs, ascending=False)

    table_rows = []
    display_cols = [
        "variant_key", "gene_name", "output_type", "biosample_name",
        "raw_score", "quantile_score",
    ]
    for _, row in sig_df.head(200).iterrows():
        cells = []
        for col in display_cols:
            val = row.get(col, "")
            if isinstance(val, float):
                val = f"{val:.6f}"
            cells.append(f"<td>{val}</td>")
        table_rows.append(f"<tr>{''.join(cells)}</tr>")

    table_html = "\n".join(table_rows) if table_rows else "<tr><td colspan='6'>No significant results</td></tr>"

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
    * {{ margin: 0; padding: 0; box-sizing: border-box; }}
    body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
           background: #f5f5f5; color: #333; padding: 2rem; }}
    .container {{ max-width: 1200px; margin: 0 auto; }}
    h1 {{ color: #1a1a2e; margin-bottom: 0.5rem; }}
    .subtitle {{ color: #666; margin-bottom: 2rem; }}
    .summary {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 1rem; margin-bottom: 2rem; }}
    .card {{ background: white; border-radius: 8px; padding: 1.5rem;
             box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}
    .card .label {{ font-size: 0.85rem; color: #666; text-transform: uppercase;
                    letter-spacing: 0.05em; }}
    .card .value {{ font-size: 2rem; font-weight: 700; color: #1a1a2e; margin-top: 0.25rem; }}
    .card .value.highlight {{ color: #e74c3c; }}
    table {{ width: 100%; border-collapse: collapse; background: white;
             border-radius: 8px; overflow: hidden; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}
    th {{ background: #1a1a2e; color: white; padding: 0.75rem 1rem; text-align: left;
          font-size: 0.85rem; text-transform: uppercase; letter-spacing: 0.05em; }}
    td {{ padding: 0.6rem 1rem; border-bottom: 1px solid #eee; font-size: 0.9rem; }}
    tr:nth-child(even) {{ background: #fafafa; }}
    tr:hover {{ background: #f0f0f0; }}
    .modalities {{ margin-bottom: 2rem; }}
    .mod-tag {{ display: inline-block; background: #e8f4fd; color: #0277bd; padding: 0.25rem 0.75rem;
                border-radius: 4px; font-size: 0.8rem; margin: 0.25rem; }}
    footer {{ margin-top: 2rem; text-align: center; color: #999; font-size: 0.8rem; }}
</style>
</head>
<body>
<div class="container">
    <h1>{title}</h1>
    <p class="subtitle">Generated by AlphaVX — AlphaGenome Variant Effect Interpreter</p>

    <div class="summary">
        <div class="card">
            <div class="label">Total Variants</div>
            <div class="value">{total_variants}</div>
        </div>
        <div class="card">
            <div class="label">Significant Hits</div>
            <div class="value highlight">{sig_count}</div>
        </div>
        <div class="card">
            <div class="label">Modalities Scored</div>
            <div class="value">{len(modalities)}</div>
        </div>
    </div>

    <div class="modalities">
        <h3>Modalities</h3>
        {''.join(f'<span class="mod-tag">{m}</span>' for m in modalities)}
    </div>

    <h3 style="margin-bottom: 1rem;">Significant Hits (quantile &gt; {quantile_threshold})</h3>
    <table>
        <thead>
            <tr>
                <th>Variant</th>
                <th>Gene</th>
                <th>Modality</th>
                <th>Tissue</th>
                <th>Raw Score</th>
                <th>Quantile Score</th>
            </tr>
        </thead>
        <tbody>
            {table_html}
        </tbody>
    </table>

    <footer>
        <p>AlphaVX &mdash; Built on AlphaGenome by Google DeepMind</p>
    </footer>
</div>
</body>
</html>"""

    report_path = output_dir / "report.html"
    report_path.write_text(html, encoding="utf-8")
    logger.info("HTML report: %s", report_path)
    return report_path


def generate_report(
    results_dir: Path,
    quantile_threshold: float = 0.995,
) -> None:
    """Read scores.csv from a results directory and generate an HTML report.

    Args:
        results_dir: Directory containing scores.csv.
        quantile_threshold: Cutoff for significance highlighting.

    Raises:
        FileNotFoundError: If scores.csv doesn't exist in results_dir.
    """
    scores_path = Path(results_dir) / "scores.csv"
    if not scores_path.exists():
        raise FileNotFoundError(f"No scores.csv found in {results_dir}")

    df = pd.read_csv(scores_path)
    generate_html_report(df, results_dir, quantile_threshold=quantile_threshold)
