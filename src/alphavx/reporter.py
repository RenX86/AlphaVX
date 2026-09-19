"""Report generation for AlphaVX — CSV, annotated VCF, and interactive HTML output."""

from __future__ import annotations

import html
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


def generate_vcf_report(
    vcf_path: Path,
    df: pd.DataFrame,
    output_dir: Path,
    quantile_threshold: float = 0.995,
) -> Path:
    """Generate an annotated VCF with AlphaVX scores in the INFO column.

    Args:
        vcf_path: Path to the original input VCF file.
        df: Complete scoring results DataFrame.
        output_dir: Directory to write the annotated VCF.
        quantile_threshold: Cutoff for significance flagging.

    Returns:
        Path to the annotated annotated.vcf file.
    """
    import gzip

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / "annotated.vcf"

    # Pre-compute variant annotations from DF
    annotations = {}
    if not df.empty and "variant_key" in df.columns:
        if "quantile_score" in df.columns:
            df["_is_sig"] = df["quantile_score"].abs() > quantile_threshold
        else:
            df["_is_sig"] = False

        for vk, group in df.groupby("variant_key"):
            max_row = (
                group.loc[group["raw_score"].abs().idxmax()]
                if "raw_score" in group.columns
                else None
            )
            max_score = max_row["raw_score"] if max_row is not None else 0.0

            sig_mods = (
                group[group["_is_sig"]]["output_type"].unique()
                if "output_type" in group.columns
                else []
            )

            annotations[vk] = {
                "max_score": max_score,
                "is_sig": len(sig_mods) > 0,
                "sig_mods": ",".join(sig_mods) if len(sig_mods) > 0 else None,
            }

    # Write annotated VCF
    opener = gzip.open(vcf_path, "rt") if str(vcf_path).endswith(".gz") else open(vcf_path)
    with opener as f_in, open(out_path, "w") as f_out:
        for line in f_in:
            line = line.strip()
            if line.startswith("##"):
                f_out.write(line + "\n")
                continue

            if line.startswith("#CHROM"):
                # Inject new INFO headers before the CHROM line
                _sig = (
                    '##INFO=<ID=AVX_SIG,Number=0,Type=Flag,Description="AlphaVX significant hit">'
                )
                _mod = (
                    "##INFO=<ID=AVX_MOD,Number=.,"
                    "Type=String,"
                    'Description="AlphaVX significant modalities">'
                )
                _max = (
                    "##INFO=<ID=AVX_MAX,Number=1,"
                    "Type=Float,"
                    'Description="AlphaVX max absolute raw score">'
                )
                f_out.write(_sig + "\n")
                f_out.write(_mod + "\n")
                f_out.write(_max + "\n")
                f_out.write(line + "\n")
                continue

            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 8:
                f_out.write(line + "\n")
                continue

            chrom = parts[0] if parts[0].startswith("chr") else f"chr{parts[0]}"
            pos = parts[1]
            ref = parts[3]
            alt = parts[4]

            # Note: For multi-allelic sites, we check the first alt allele for simplicity
            first_alt = alt.split(",")[0]
            vk = f"{chrom}:{pos}:{ref}>{first_alt}"

            if vk in annotations:
                ann = annotations[vk]
                info = parts[7]
                new_info = []
                if info != ".":
                    new_info.append(info)

                new_info.append(f"AVX_MAX={ann['max_score']:.5f}")
                if ann["is_sig"]:
                    new_info.append("AVX_SIG")
                    if ann["sig_mods"]:
                        new_info.append(f"AVX_MOD={ann['sig_mods']}")

                parts[7] = ";".join(new_info)

            f_out.write("\t".join(parts) + "\n")

    logger.info("Annotated VCF: %s", out_path)
    return out_path


def generate_html_report(
    df: pd.DataFrame,
    output_dir: Path,
    title: str = "AlphaVX Variant Report",
    quantile_threshold: float = 0.995,
) -> Path:
    """Generate an interactive HTML report from scoring results.

    Uses DataTables.js for searchable/sortable/paginated tables and
    Plotly.js for interactive heatmap, volcano plot, and modality
    breakdown charts.  The report is a single self-contained HTML file;
    JS libraries are loaded from CDN at view-time.

    Args:
        df: Complete scoring results DataFrame.
        output_dir: Directory to write the report.
        title: Report title.
        quantile_threshold: Cutoff for significance highlighting.

    Returns:
        Path to the generated report.html file.
    """
    import json as _json
    from datetime import datetime, timezone

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Compute summary stats
    # ------------------------------------------------------------------
    if "quantile_score" in df.columns:
        sig_mask = df["quantile_score"].abs() > quantile_threshold
    else:
        sig_mask = pd.Series([False] * len(df))

    total_variants = df["variant_key"].nunique() if "variant_key" in df.columns else 0
    sig_variant_count = (
        df.loc[sig_mask, "variant_key"].nunique() if "variant_key" in df.columns else 0
    )
    total_scores = len(df)
    sig_scores = int(sig_mask.sum())
    modalities = sorted(df["output_type"].unique().tolist()) if "output_type" in df.columns else []

    # ------------------------------------------------------------------
    # Build JSON data for DataTables + Plotly
    # ------------------------------------------------------------------
    sig_df = df[sig_mask].copy() if sig_mask.any() else pd.DataFrame()
    if not sig_df.empty and "quantile_score" in sig_df.columns:
        sig_df = sig_df.sort_values("quantile_score", key=abs, ascending=False)

    display_cols = [
        "variant_key",
        "gene_name",
        "output_type",
        "biosample_name",
        "raw_score",
        "quantile_score",
    ]
    table_data = []
    for _, row in sig_df.iterrows():
        record = []
        for col in display_cols:
            val = row.get(col, "")
            if isinstance(val, float):
                val = round(val, 6)
            else:
                val = html.escape(str(val))
            record.append(val)
        table_data.append(record)

    table_json = _json.dumps(table_data)

    # Heatmap data: pivot variant x modality -> max abs raw_score
    heatmap_json = "{}"
    if (
        not sig_df.empty
        and "variant_key" in sig_df.columns
        and "output_type" in sig_df.columns
        and "raw_score" in sig_df.columns
    ):
        pivot = sig_df.pivot_table(
            index="variant_key",
            columns="output_type",
            values="raw_score",
            aggfunc=lambda x: x.iloc[x.abs().argmax()],
        )
        heatmap_json = _json.dumps(
            {
                "variants": pivot.index.tolist(),
                "modalities": pivot.columns.tolist(),
                "z": [[round(v, 4) if pd.notna(v) else None for v in row] for row in pivot.values],
            }
        )

    # Volcano data: raw_score vs -log10(1 - |quantile|)
    volcano_json = "[]"
    if not sig_df.empty and "raw_score" in sig_df.columns and "quantile_score" in sig_df.columns:
        import math

        volcano_records = []
        for _, row in sig_df.iterrows():
            q = abs(row.get("quantile_score", 0))
            neg_log = -math.log10(max(1 - q, 1e-10))
            volcano_records.append(
                {
                    "x": round(row.get("raw_score", 0), 6),
                    "y": round(neg_log, 4),
                    "variant": html.escape(str(row.get("variant_key", ""))),
                    "gene": html.escape(str(row.get("gene_name", ""))),
                    "modality": html.escape(str(row.get("output_type", ""))),
                }
            )
        volcano_json = _json.dumps(volcano_records)

    # Modality breakdown: count of sig hits per variant per modality
    breakdown_json = "{}"
    if not sig_df.empty and "variant_key" in sig_df.columns and "output_type" in sig_df.columns:
        bd = sig_df.groupby(["variant_key", "output_type"]).size().unstack(fill_value=0)
        breakdown_json = _json.dumps(
            {
                "variants": bd.index.tolist(),
                "modalities": bd.columns.tolist(),
                "counts": bd.values.tolist(),
            }
        )

    # ------------------------------------------------------------------
    # Escape user-supplied strings
    # ------------------------------------------------------------------
    escaped_title = html.escape(title)
    escaped_mod_tags = "".join(f'<span class="mod-tag">{html.escape(m)}</span>' for m in modalities)
    generated_at = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    # ------------------------------------------------------------------
    # HTML template
    # ------------------------------------------------------------------
    # Note: double braces {{ }} are Python f-string escapes for literal { }
    html_out = _build_html_template(
        escaped_title=escaped_title,
        total_variants=total_variants,
        sig_variant_count=sig_variant_count,
        total_scores=total_scores,
        sig_scores=sig_scores,
        modalities=modalities,
        escaped_mod_tags=escaped_mod_tags,
        table_json=table_json,
        heatmap_json=heatmap_json,
        volcano_json=volcano_json,
        breakdown_json=breakdown_json,
        quantile_threshold=quantile_threshold,
        generated_at=generated_at,
    )

    report_path = output_dir / "report.html"
    report_path.write_text(html_out, encoding="utf-8")
    logger.info("HTML report: %s", report_path)
    return report_path


def _build_html_template(
    *,
    escaped_title: str,
    total_variants: int,
    sig_variant_count: int,
    total_scores: int,
    sig_scores: int,
    modalities: list[str],
    escaped_mod_tags: str,
    table_json: str,
    heatmap_json: str,
    volcano_json: str,
    breakdown_json: str,
    quantile_threshold: float,
    generated_at: str,
) -> str:
    """Build the full interactive HTML string.

    Separated from generate_html_report for readability.  Uses Python
    string formatting with {{ }} for literal braces needed by CSS/JS.
    """
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escaped_title}</title>

<!-- jQuery + DataTables -->
<link rel="stylesheet" href="https://cdn.datatables.net/1.13.7/css/jquery.dataTables.min.css">
<link rel="stylesheet" href="https://cdn.datatables.net/buttons/2.4.2/css/buttons.dataTables.min.css">

<!-- Plotly -->
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js" charset="utf-8"></script>

<style>
    :root {{
        --bg: #f5f5f5; --bg-card: #fff; --text: #333; --text-muted: #666;
        --header-bg: #1a1a2e; --header-text: #fff;
        --accent: #e74c3c; --accent2: #0277bd;
        --border: #eee; --hover: #f0f0f0;
    }}
    [data-theme="dark"] {{
        --bg: #121212; --bg-card: #1e1e1e; --text: #e0e0e0; --text-muted: #999;
        --header-bg: #0d0d2b; --header-text: #e0e0e0;
        --border: #333; --hover: #2a2a2a;
    }}
    * {{ margin: 0; padding: 0; box-sizing: border-box; }}
    body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
           background: var(--bg); color: var(--text); padding: 2rem;
           transition: background 0.3s, color 0.3s; }}
    .container {{ max-width: 1400px; margin: 0 auto; }}

    /* Header */
    .header {{ display: flex; justify-content: space-between; align-items: center;
               margin-bottom: 0.5rem; }}
    h1 {{ color: var(--text); font-size: 1.8rem; }}
    .subtitle {{ color: var(--text-muted); margin-bottom: 2rem; }}
    .theme-toggle {{ background: var(--bg-card); border: 1px solid var(--border);
                     border-radius: 6px; padding: 0.4rem 0.8rem; cursor: pointer;
                     font-size: 0.85rem; color: var(--text); }}

    /* Summary cards */
    .summary {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
                gap: 1rem; margin-bottom: 2rem; }}
    .card {{ background: var(--bg-card); border-radius: 8px; padding: 1.5rem;
             box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}
    .card .label {{ font-size: 0.8rem; color: var(--text-muted); text-transform: uppercase;
                    letter-spacing: 0.05em; }}
    .card .value {{ font-size: 1.8rem; font-weight: 700; color: var(--text);
                    margin-top: 0.25rem; }}
    .card .value.highlight {{ color: var(--accent); }}
    .modalities {{ margin-bottom: 2rem; }}
    .mod-tag {{ display: inline-block; background: #e8f4fd; color: var(--accent2);
                padding: 0.25rem 0.75rem; border-radius: 4px; font-size: 0.8rem;
                margin: 0.25rem; }}

    /* Tabs */
    .tabs {{ display: flex; gap: 0; margin-bottom: 0;
             border-bottom: 2px solid var(--border); }}
    .tab-btn {{ padding: 0.75rem 1.5rem; cursor: pointer; border: none; background: none;
                font-size: 0.95rem; color: var(--text-muted);
                border-bottom: 2px solid transparent; margin-bottom: -2px;
                transition: all 0.2s; }}
    .tab-btn:hover {{ color: var(--text); }}
    .tab-btn.active {{ color: var(--accent2); border-bottom-color: var(--accent2);
                       font-weight: 600; }}
    .tab-content {{ display: none; padding: 1.5rem 0; }}
    .tab-content.active {{ display: block; }}

    /* Charts */
    .chart-container {{ background: var(--bg-card); border-radius: 8px; padding: 1rem;
                        box-shadow: 0 1px 3px rgba(0,0,0,0.1); margin-bottom: 1.5rem; }}
    .charts-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 1.5rem; }}

    /* DataTables overrides */
    .dataTables_wrapper {{ background: var(--bg-card); border-radius: 8px; padding: 1rem;
                           box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}
    table.dataTable thead th {{ background: var(--header-bg) !important;
                                color: var(--header-text) !important; }}
    table.dataTable tbody td {{ border-bottom: 1px solid var(--border); }}
    table.dataTable tbody tr:hover {{ background: var(--hover) !important; }}

    footer {{ margin-top: 2rem; text-align: center; color: var(--text-muted);
              font-size: 0.8rem; }}

    @media (max-width: 900px) {{
        .charts-grid {{ grid-template-columns: 1fr; }}
    }}
</style>
</head>
<body>
<div class="container">

    <div class="header">
        <h1>{escaped_title}</h1>
        <button class="theme-toggle" onclick="toggleTheme()">&#x1F319; Dark Mode</button>
    </div>
    <p class="subtitle">Generated by AlphaVX &mdash; AlphaGenome Variant Effect Interpreter</p>

    <!-- Summary Cards -->
    <div class="summary">
        <div class="card">
            <div class="label">Total Variants</div>
            <div class="value">{total_variants}</div>
        </div>
        <div class="card">
            <div class="label">Significant Variants</div>
            <div class="value highlight">{sig_variant_count}</div>
        </div>
        <div class="card">
            <div class="label">Total Scores</div>
            <div class="value">{total_scores:,}</div>
        </div>
        <div class="card">
            <div class="label">Significant Scores</div>
            <div class="value highlight">{sig_scores:,}</div>
        </div>
        <div class="card">
            <div class="label">Modalities</div>
            <div class="value">{len(modalities)}</div>
        </div>
    </div>

    <div class="modalities">
        <h3>Modalities Scored</h3>
        {escaped_mod_tags}
    </div>

    <!-- Tabs -->
    <div class="tabs" role="tablist">
        <button class="tab-btn active" onclick="switchTab('overview')" role="tab">Overview</button>
        <button class="tab-btn" onclick="switchTab('table')" role="tab">Data Table</button>
        <button class="tab-btn" onclick="switchTab('charts')" role="tab">Charts</button>
    </div>

    <!-- Tab: Overview (heatmap) -->
    <div id="tab-overview" class="tab-content active" role="tabpanel">
        <div class="chart-container">
            <div id="heatmap" style="width:100%; min-height:400px;"></div>
        </div>
    </div>

    <!-- Tab: Data Table -->
    <div id="tab-table" class="tab-content" role="tabpanel">
        <table id="sig-table" class="display" style="width:100%">
            <thead>
                <tr>
                    <th>Variant</th><th>Gene</th><th>Modality</th>
                    <th>Tissue</th><th>Raw Score</th><th>Quantile Score</th>
                </tr>
            </thead>
        </table>
    </div>

    <!-- Tab: Charts (volcano + breakdown) -->
    <div id="tab-charts" class="tab-content" role="tabpanel">
        <div class="charts-grid">
            <div class="chart-container">
                <div id="volcano" style="width:100%; min-height:450px;"></div>
            </div>
            <div class="chart-container">
                <div id="breakdown" style="width:100%; min-height:450px;"></div>
            </div>
        </div>
    </div>

    <footer>
        <p>AlphaVX &mdash; Built on AlphaGenome by Google DeepMind &bull; {generated_at}</p>
    </footer>
</div>

<!-- jQuery + DataTables JS -->
<script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
<script src="https://cdn.datatables.net/1.13.7/js/jquery.dataTables.min.js"></script>
<script src="https://cdn.datatables.net/buttons/2.4.2/js/dataTables.buttons.min.js"></script>
<script src="https://cdn.datatables.net/buttons/2.4.2/js/buttons.html5.min.js"></script>

<script>
// ---- Embedded data ----
var TABLE_DATA = {table_json};
var HEATMAP_DATA = {heatmap_json};
var VOLCANO_DATA = {volcano_json};
var BREAKDOWN_DATA = {breakdown_json};
var THRESHOLD = {quantile_threshold};

// ---- Theme toggle ----
function toggleTheme() {{
    var el = document.documentElement;
    var btn = document.querySelector('.theme-toggle');
    if (el.getAttribute('data-theme') === 'dark') {{
        el.removeAttribute('data-theme');
        btn.innerHTML = '&#x1F319; Dark Mode';
    }} else {{
        el.setAttribute('data-theme', 'dark');
        btn.innerHTML = '&#x2600;&#xFE0F; Light Mode';
    }}
    renderCharts();
}}

// ---- Tab switching ----
function switchTab(name) {{
    document.querySelectorAll('.tab-content').forEach(function(el) {{
        el.classList.remove('active');
    }});
    document.querySelectorAll('.tab-btn').forEach(function(el) {{
        el.classList.remove('active');
    }});
    document.getElementById('tab-' + name).classList.add('active');
    var labels = {{ overview: 'Overview', table: 'Data Table', charts: 'Charts' }};
    document.querySelectorAll('.tab-btn').forEach(function(el) {{
        if (el.textContent.trim() === labels[name]) el.classList.add('active');
    }});
    setTimeout(function() {{ window.dispatchEvent(new Event('resize')); }}, 100);
}}

// ---- DataTables init ----
$(document).ready(function() {{
    $('#sig-table').DataTable({{
        data: TABLE_DATA,
        pageLength: 25,
        lengthMenu: [10, 25, 50, 100],
        order: [[5, 'desc']],
        dom: 'Blfrtip',
        buttons: ['csv', 'copy'],
        columns: [
            {{ title: 'Variant' }},
            {{ title: 'Gene' }},
            {{ title: 'Modality' }},
            {{ title: 'Tissue' }},
            {{ title: 'Raw Score',
              render: function(d) {{ return typeof d === 'number' ? d.toFixed(6) : d; }} }},
            {{ title: 'Quantile',
              render: function(d) {{ return typeof d === 'number' ? d.toFixed(6) : d; }} }},
        ],
    }});
}});

// ---- Plotly charts ----
function getPlotlyTheme() {{
    var dark = document.documentElement.getAttribute('data-theme') === 'dark';
    return {{
        paper_bgcolor: dark ? '#1e1e1e' : '#fff',
        plot_bgcolor:  dark ? '#1e1e1e' : '#fff',
        font: {{ color: dark ? '#e0e0e0' : '#333' }},
    }};
}}

function renderCharts() {{
    var t = getPlotlyTheme();

    // -- Heatmap --
    if (HEATMAP_DATA.z && HEATMAP_DATA.z.length > 0) {{
        Plotly.newPlot('heatmap', [{{
            z: HEATMAP_DATA.z,
            x: HEATMAP_DATA.modalities,
            y: HEATMAP_DATA.variants,
            type: 'heatmap',
            colorscale: 'RdBu',
            zmid: 0,
            colorbar: {{ title: 'Raw Score' }},
            hovertemplate:
                'Variant: %{{y}}<br>Modality: %{{x}}<br>Score: %{{z:.4f}}<extra></extra>',
        }}], Object.assign({{
            title: 'Significant Variant Effects by Modality',
            xaxis: {{ title: 'Modality', tickangle: -45 }},
            yaxis: {{ title: 'Variant', automargin: true }},
            margin: {{ l: 180, b: 120, t: 50 }},
        }}, t), {{ responsive: true }});
    }} else {{
        document.getElementById('heatmap').innerHTML =
            '<p style="text-align:center;padding:3rem;color:#999;">'
            + 'No significant results to display</p>';
    }}

    // -- Volcano plot --
    if (VOLCANO_DATA.length > 0) {{
        Plotly.newPlot('volcano', [{{
            x: VOLCANO_DATA.map(function(d) {{ return d.x; }}),
            y: VOLCANO_DATA.map(function(d) {{ return d.y; }}),
            text: VOLCANO_DATA.map(function(d) {{ return d.gene + ' (' + d.modality + ')'; }}),
            customdata: VOLCANO_DATA.map(function(d) {{ return d.variant; }}),
            mode: 'markers',
            type: 'scatter',
            marker: {{
                size: 5,
                color: VOLCANO_DATA.map(function(d) {{ return d.x; }}),
                colorscale: 'RdBu', cmid: 0, opacity: 0.7,
            }},
            hovertemplate:
                '%{{text}}<br>Variant: %{{customdata}}<br>'
                + 'Raw: %{{x:.4f}}<br>-log10(1-|q|): %{{y:.2f}}<extra></extra>',
        }}], Object.assign({{
            title: 'Volcano Plot — Effect Size vs Significance',
            xaxis: {{ title: 'Raw Score (log-fold change)', zeroline: true }},
            yaxis: {{ title: '-log\\u2081\\u2080(1 - |quantile|)' }},
            margin: {{ t: 50 }},
        }}, t), {{ responsive: true }});
    }} else {{
        document.getElementById('volcano').innerHTML =
            '<p style="text-align:center;padding:3rem;color:#999;">No data</p>';
    }}

    // -- Modality breakdown --
    if (BREAKDOWN_DATA.counts && BREAKDOWN_DATA.counts.length > 0) {{
        var traces = [];
        for (var m = 0; m < BREAKDOWN_DATA.modalities.length; m++) {{
            traces.push({{
                x: BREAKDOWN_DATA.variants,
                y: BREAKDOWN_DATA.counts.map(function(row) {{ return row[m]; }}),
                name: BREAKDOWN_DATA.modalities[m],
                type: 'bar',
            }});
        }}
        Plotly.newPlot('breakdown', traces, Object.assign({{
            title: 'Significant Hits per Variant by Modality',
            barmode: 'stack',
            xaxis: {{ title: 'Variant', tickangle: -45 }},
            yaxis: {{ title: 'Count of Significant Scores' }},
            margin: {{ b: 150, t: 50 }},
            legend: {{ orientation: 'h', y: -0.35 }},
        }}, t), {{ responsive: true }});
    }} else {{
        document.getElementById('breakdown').innerHTML =
            '<p style="text-align:center;padding:3rem;color:#999;">No data</p>';
    }}
}}

// Render on load
renderCharts();
</script>
</body>
</html>"""


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
