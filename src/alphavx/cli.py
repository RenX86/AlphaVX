"""AlphaVX command-line interface.

Usage:
    alphavx score input.vcf -o results/
    alphavx query chr17:7674220:G>A
    alphavx report results/
    alphavx cache stats
    alphavx cache clear
"""

from __future__ import annotations

import logging
from pathlib import Path

import typer

from . import __version__
from .cache import ResultCache
from .config import load_config
from .reporter import (
    generate_csv_report,
    generate_html_report,
    generate_report,
    generate_vcf_report,
)
from .scorer import VariantScorer
from .vcf_parser import parse_variant_string, parse_vcf

app = typer.Typer(
    name="alphavx",
    help="AlphaVX — AlphaGenome Variant Effect Interpreter.\n\n"
    "Batch-score genetic variants against AlphaGenome and generate "
    "interpretable multi-modal effect reports.",
    add_completion=False,
)

cache_app = typer.Typer(help="Manage the variant scoring cache.")
app.add_typer(cache_app, name="cache")

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-7s | %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"alphavx {__version__}")
        raise typer.Exit()


def _get_console():
    """Get a Rich console for pretty output, or fall back to plain print."""
    try:
        from rich.console import Console

        return Console()
    except ImportError:
        return None


@app.callback()
def main(
    version: bool | None = typer.Option(
        None,
        "--version",
        "-v",
        callback=_version_callback,
        is_eager=True,
        help="Show version and exit.",
    ),
) -> None:
    \"\"\"AlphaVX — AlphaGenome Variant Effect Interpreter.\"\"\"


@app.command()
def configure() -> None:
    \"\"\"Interactively configure AlphaVX API keys globally.\"\"\"
    typer.echo("AlphaVX Configuration")
    typer.echo("---------------------")
    api_key = typer.prompt("Enter your AlphaGenome API Key", hide_input=True)

    config_dir = Path.home() / ".alphavx"
    config_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        env_file = config_dir / ".env"
        with open(env_file, "w") as f:
            f.write(f"ALPHAVX_API_KEY={api_key}\\n")
            
        typer.secho(f"\\n✅ API Key securely saved to {env_file}", fg=typer.colors.GREEN)
        typer.echo("You can now run 'alphavx score' from any directory.")
        
    except Exception as e:
        typer.secho(f"Failed to save configuration: {e}", fg=typer.colors.RED, err=True)


@app.command()
def init() -> None:
    \"\"\"Create a sample VCF file in the current directory for testing.\"\"\"
    sample_vcf = Path.cwd() / "clinvar_example.vcf"
    
    # A tiny VCF with a few known pathogenic variants
    vcf_content = \"\"\"##fileformat=VCFv4.2
##source=AlphaVX_Example
##contig=<ID=chr17,length=83257441>
##contig=<ID=chr7,length=159345973>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr17\t43094895\t.\tC\tT\t.\t.\tGENE=BRCA1;CLIN_SIG=Pathogenic
chr7\t117559590\t.\tT\tC\t.\t.\tGENE=CFTR;CLIN_SIG=Pathogenic
chr17\t7674220\t.\tG\tA\t.\t.\tGENE=TP53;CLIN_SIG=Pathogenic
\"\"\"
    
    try:
        with open(sample_vcf, "w") as f:
            f.write(vcf_content)
        typer.secho(f"✅ Created {sample_vcf.name}", fg=typer.colors.GREEN)
        typer.echo("\\nTry running your first batch score:")
        typer.secho(f"  alphavx score {sample_vcf.name} -o results/", fg=typer.colors.CYAN)
    except Exception as e:
        typer.secho(f"Failed to create example: {e}", fg=typer.colors.RED, err=True)


@app.command()
def score(
    vcf_path: Path = typer.Argument(..., help="Path to VCF file containing variants to score."),
    output: Path = typer.Option("results", "--output", "-o", help="Output directory for results."),
    genes: str | None = typer.Option(
        None,
        "--genes",
        "-g",
        help="Comma-separated gene list to filter results (e.g., BRCA1,TP53,CFTR).",
    ),
    config_path: Path | None = typer.Option(
        None,
        "--config",
        "-c",
        help="Path to alphavx.yaml configuration file.",
    ),
    no_cache: bool = typer.Option(False, "--no-cache", help="Disable result caching."),
) -> None:
    """Score all variants in a VCF file against AlphaGenome."""
    _get_console()

    # Load config
    try:
        config = load_config(config_path)
    except ValueError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1)

    config.output_dir = output
    if no_cache:
        config.cache_enabled = False

    # Parse VCF
    try:
        records = parse_vcf(vcf_path)
    except (FileNotFoundError, ValueError) as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1)

    typer.echo(f"Parsed {len(records)} variants from {vcf_path}")

    # Setup cache
    cache = None
    if config.cache_enabled:
        cache = ResultCache(output / "cache", config_hash=config.scoring_fingerprint)
        stats = cache.stats()
        if stats["count"] > 0:
            typer.echo(f"Cache: {stats['count']} variants already cached")

    # Score
    scorer = VariantScorer(config)

    def progress(i: int, total: int, record) -> None:
        typer.echo(f"[{i + 1}/{total}] Scoring {record.key}...")

    df = scorer.score_batch(records, cache=cache, progress_callback=progress)

    if df.empty:
        typer.echo("No results — all variants failed or returned empty scores.")
        raise typer.Exit(1)

    # Filter by genes if specified
    if genes:
        gene_list = [g.strip().upper() for g in genes.split(",")]
        gene_col = "gene_name" if "gene_name" in df.columns else None
        if gene_col:
            before = len(df)
            df = df[df[gene_col].str.upper().isin(gene_list)]
            typer.echo(f"Gene filter: {before} → {len(df)} rows (genes: {', '.join(gene_list)})")

    # Generate reports
    generate_csv_report(df, output, quantile_threshold=config.quantile_threshold)
    generate_html_report(df, output, quantile_threshold=config.quantile_threshold)

    # Generate annotated VCF
    try:
        generate_vcf_report(vcf_path, df, output, quantile_threshold=config.quantile_threshold)
    except Exception as e:
        logger.warning("Annotated VCF generation failed: %s", e)

    # Generate plots
    try:
        from .plots import plot_summary_heatmap, plot_variant_detail

        plots_dir = output / "plots"
        plot_summary_heatmap(df, plots_dir / "summary_heatmap.png", config.quantile_threshold)

        # Generate per-variant detail plots
        variant_keys = df["variant_key"].unique() if "variant_key" in df.columns else []
        for vk in variant_keys:
            safe_name = str(vk).replace(":", "_").replace(">", "_")
            plot_variant_detail(
                df,
                vk,
                plots_dir / "per_variant" / f"{safe_name}.png",
                config.quantile_threshold,
            )
        if len(variant_keys):
            logger.info("Generated %d per-variant plots", len(variant_keys))
    except Exception as e:
        logger.warning("Plot generation failed: %s", e)

    # Print summary
    sig_count = 0
    if "quantile_score" in df.columns:
        sig_count = (df["quantile_score"].abs() > config.quantile_threshold).sum()

    typer.echo("")
    typer.echo("═" * 50)
    n_scored = df["variant_key"].nunique() if "variant_key" in df.columns else "?"
    typer.echo(f"  Variants scored:    {n_scored}")
    typer.echo(f"  Significant hits:   {sig_count}")
    typer.echo(f"  Results:            {output / 'scores.csv'}")
    typer.echo(f"  Report:             {output / 'report.html'}")
    typer.echo("═" * 50)


@app.command()
def query(
    variant: str = typer.Argument(
        ...,
        help="Variant to score (e.g., chr17:7674220:G>A or chr17:7674220:G:A).",
    ),
    output: Path | None = typer.Option(
        None,
        "--output",
        "-o",
        help="Optional directory to save full CSV results.",
    ),
    config_path: Path | None = typer.Option(
        None,
        "--config",
        "-c",
        help="Path to alphavx.yaml configuration file.",
    ),
) -> None:
    """Score a single variant and display results."""
    # Parse variant
    try:
        record = parse_variant_string(variant)
    except ValueError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1)

    # Load config
    try:
        config = load_config(config_path)
    except ValueError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1)

    typer.echo(f"Scoring {record.key}...")

    # Score
    scorer = VariantScorer(config)
    df = scorer.score_variant(record)

    if df.empty:
        typer.echo("No results returned for this variant.")
        raise typer.Exit(1)

    # Display significant results
    if "quantile_score" in df.columns:
        sig = df[df["quantile_score"].abs() > config.quantile_threshold]
    else:
        sig = df.head(0)

    typer.echo(f"\nTotal scores: {len(df)}")
    typer.echo(f"Significant:  {len(sig)}")

    if not sig.empty:
        header = (
            f"{'Variant':<25} {'Gene':<12} {'Modality':<20}"
            f" {'Tissue':<25} {'Raw':<12} {'Quantile':<10}"
        )
        typer.echo(f"\n{header}")
        typer.echo("─" * 104)
        for _, row in (
            sig.sort_values("quantile_score", key=abs, ascending=False).head(30).iterrows()
        ):
            vk = str(row.get("variant_key", ""))[:24]
            gene = str(row.get("gene_name", ""))[:11]
            mod = str(row.get("output_type", ""))[:19]
            tissue = str(row.get("biosample_name", ""))[:24]
            raw = f"{row.get('raw_score', 0):.6f}"
            quant = f"{row.get('quantile_score', 0):.6f}"
            typer.echo(f"{vk:<25} {gene:<12} {mod:<20} {tissue:<25} {raw:<12} {quant:<10}")

    # Save if output specified
    if output:
        generate_csv_report(df, output, quantile_threshold=config.quantile_threshold)
        typer.echo(f"\nFull results saved to {output / 'scores.csv'}")


@app.command()
def report(
    results_dir: Path = typer.Argument(
        ...,
        help="Directory containing scores.csv to generate report from.",
    ),
) -> None:
    """Generate an HTML report from existing scoring results."""
    try:
        generate_report(results_dir)
        typer.echo(f"Report generated: {results_dir / 'report.html'}")
    except FileNotFoundError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(1)


@cache_app.command("stats")
def cache_stats(
    cache_dir: Path = typer.Option(
        "results/cache",
        "--dir",
        "-d",
        help="Cache directory.",
    ),
) -> None:
    """Show cache statistics."""
    cache = ResultCache(cache_dir)
    stats = cache.stats()
    typer.echo(f"Cached variants: {stats['count']}")
    typer.echo(f"Cache size:      {stats['cache_size_bytes'] / 1024:.1f} KB")
    typer.echo(f"Cache location:  {cache.db_path}")


@cache_app.command("clear")
def cache_clear(
    cache_dir: Path = typer.Option(
        "results/cache",
        "--dir",
        "-d",
        help="Cache directory.",
    ),
    force: bool = typer.Option(False, "--force", "-f", help="Skip confirmation prompt."),
) -> None:
    """Clear all cached results."""
    cache = ResultCache(cache_dir)
    stats = cache.stats()

    if stats["count"] == 0:
        typer.echo("Cache is already empty.")
        return

    if not force:
        confirm = typer.confirm(f"Delete {stats['count']} cached variant scores?")
        if not confirm:
            typer.echo("Cancelled.")
            return

    cache.clear()
    typer.echo(f"Cleared {stats['count']} cached entries.")
