"""Visualization module for AlphaVX — summary heatmaps and per-variant plots."""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # Headless rendering — no display required

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

logger = logging.getLogger(__name__)


def plot_summary_heatmap(
    df: pd.DataFrame,
    output_path: Path,
    quantile_threshold: float = 0.995,
) -> Path:
    """Create a heatmap of significant variants × modalities.

    Rows are variants, columns are modalities (output_type), values
    are raw_score. Only significant variants (|quantile_score| > threshold)
    are included.

    Args:
        df: Full scoring results DataFrame.
        output_path: Path to save the PNG file.
        quantile_threshold: Cutoff for including variants.

    Returns:
        Path to the saved figure.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Filter to significant
    if "quantile_score" not in df.columns or "output_type" not in df.columns:
        _save_empty_figure(output_path, "No data available for heatmap")
        return output_path

    sig = df[df["quantile_score"].abs() > quantile_threshold].copy()

    if sig.empty:
        _save_empty_figure(output_path, "No significant variants found")
        return output_path

    # Build pivot table: variant × modality → raw_score (take max abs per cell)
    if "variant_key" not in sig.columns:
        sig["variant_key"] = sig.apply(
            lambda r: (
                f"{r.get('chrom', '')}:{r.get('pos', '')}:"
                f"{r.get('ref', '')}>{r.get('alt', '')}"
            ),
            axis=1,
        )

    pivot = sig.pivot_table(
        index="variant_key",
        columns="output_type",
        values="raw_score",
        aggfunc=lambda x: x.iloc[x.abs().argmax()],  # max absolute score
    )

    # Plot
    fig_height = max(4, len(pivot) * 0.5 + 2)
    fig_width = max(8, len(pivot.columns) * 1.2 + 3)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    sns.heatmap(
        pivot,
        cmap="RdBu_r",
        center=0,
        annot=True,
        fmt=".3f",
        linewidths=0.5,
        ax=ax,
        cbar_kws={"label": "Raw Score (log-fold change)"},
    )
    ax.set_title("AlphaVX — Significant Variant Effects by Modality", fontsize=14, pad=15)
    ax.set_xlabel("Modality", fontsize=11)
    ax.set_ylabel("Variant", fontsize=11)
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    plt.tight_layout()

    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("Summary heatmap: %s", output_path)
    return output_path


def plot_variant_detail(
    df: pd.DataFrame,
    variant_key: str,
    output_path: Path,
    quantile_threshold: float = 0.995,
) -> Path:
    """Create a bar chart of scores across modalities for a single variant.

    Bars are colored by significance (red if significant, grey otherwise).

    Args:
        df: Full scoring results DataFrame.
        variant_key: Variant identifier (chr:pos:ref>alt).
        output_path: Path to save the PNG file.
        quantile_threshold: Cutoff for significance coloring.

    Returns:
        Path to the saved figure.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    vdf = df[df["variant_key"] == variant_key].copy()

    if vdf.empty:
        _save_empty_figure(output_path, f"No data for variant {variant_key}")
        return output_path

    # Aggregate: max absolute raw_score per modality
    if "output_type" not in vdf.columns or "raw_score" not in vdf.columns:
        _save_empty_figure(output_path, f"Missing columns for {variant_key}")
        return output_path

    agg = (
        vdf.groupby("output_type")
        .agg(
            raw_score=("raw_score", lambda x: x.iloc[x.abs().argmax()]),
            quantile_score=("quantile_score", lambda x: x.iloc[x.abs().argmax()]),
        )
        .reset_index()
    )
    agg = agg.sort_values("raw_score", key=abs, ascending=True)

    # Color by significance
    colors = [
        "#e74c3c" if abs(q) > quantile_threshold else "#95a5a6" for q in agg["quantile_score"]
    ]

    fig, ax = plt.subplots(figsize=(10, max(4, len(agg) * 0.5 + 1)))
    ax.barh(agg["output_type"], agg["raw_score"], color=colors, edgecolor="white")

    ax.axvline(x=0, color="black", linewidth=0.8)
    ax.set_xlabel("Raw Score (log-fold change)", fontsize=11)
    ax.set_title(f"AlphaVX — {variant_key}", fontsize=13, pad=10)

    # Add significance legend
    from matplotlib.patches import Patch

    legend_elements = [
        Patch(facecolor="#e74c3c", label=f"Significant (|q| > {quantile_threshold})"),
        Patch(facecolor="#95a5a6", label="Not significant"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=9)

    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("Variant detail plot: %s", output_path)
    return output_path


def _save_empty_figure(output_path: Path, message: str) -> None:
    """Save a placeholder figure with a text message."""
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.text(
        0.5,
        0.5,
        message,
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=14,
        color="#999",
    )
    ax.set_axis_off()
    fig.savefig(output_path, dpi=100, bbox_inches="tight")
    plt.close(fig)
