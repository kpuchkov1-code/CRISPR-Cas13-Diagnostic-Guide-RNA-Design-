"""Visualization functions for guide RNA analysis."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns


sns.set_theme(style="whitegrid", palette="deep", font_scale=1.1)


def plot_conservation_heatmap(
    position_scores: np.ndarray,
    candidates: pd.DataFrame = None,
    figsize: tuple = (14, 3),
) -> plt.Figure:
    """Plot conservation score across the alignment as a heatmap.

    Optionally marks candidate guide positions with vertical lines.

    Args:
        position_scores: Per-position conservation array.
        candidates: Optional DataFrame with 'position' column to mark.
        figsize: Figure dimensions.

    Returns:
        matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Reshape for heatmap display
    data = position_scores.reshape(1, -1)
    cmap = sns.color_palette("YlOrRd", as_cmap=True)

    im = ax.imshow(
        data, aspect="auto", cmap=cmap, vmin=0, vmax=1,
        extent=[0, len(position_scores), 0, 1],
    )

    if candidates is not None and not candidates.empty:
        for _, row in candidates.iterrows():
            ax.axvline(
                x=row["position"], color="#2196F3", alpha=0.6,
                linewidth=0.8, linestyle="--",
            )

    ax.set_xlabel("Alignment Position")
    ax.set_yticks([])
    ax.set_title("Sequence Conservation Across Alignment", fontweight="bold")

    cbar = plt.colorbar(im, ax=ax, orientation="horizontal", pad=0.25, shrink=0.5)
    cbar.set_label("Conservation Score")

    plt.tight_layout()
    return fig


def plot_score_breakdown(
    ranked_df: pd.DataFrame,
    top_n: int = 10,
    figsize: tuple = (12, 6),
) -> plt.Figure:
    """Stacked bar chart showing composite score breakdown for top guides.

    Args:
        ranked_df: DataFrame from rank_guides.
        top_n: Number of guides to display.
        figsize: Figure dimensions.

    Returns:
        matplotlib Figure.
    """
    df = ranked_df.head(top_n).copy()

    if df.empty:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "No candidates to display", ha="center", va="center")
        return fig

    score_cols = [
        ("score_conservation", "Conservation", "#4CAF50"),
        ("score_offtarget", "Off-Target Safety", "#2196F3"),
        ("score_gc", "GC Optimality", "#FF9800"),
        ("score_homopolymer", "Homopolymer", "#9C27B0"),
        ("score_selfcomp", "Self-Folding", "#F44336"),
    ]

    weights = [0.30, 0.25, 0.15, 0.15, 0.15]

    fig, ax = plt.subplots(figsize=figsize)

    labels = [f"Guide {r}" for r in df["rank"]]
    x = np.arange(len(labels))
    bottom = np.zeros(len(df))

    for (col, label, color), w in zip(score_cols, weights):
        if col in df.columns:
            values = df[col].values * w
            ax.bar(x, values, bottom=bottom, label=f"{label} ({int(w*100)}%)",
                   color=color, alpha=0.85, edgecolor="white", linewidth=0.5)
            bottom += values

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("Weighted Score Contribution")
    ax.set_title("Top Guide RNAs — Composite Score Breakdown", fontweight="bold")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
    ax.set_ylim(0, 1.05)

    plt.tight_layout()
    return fig


def plot_gc_distribution(
    candidates: pd.DataFrame,
    figsize: tuple = (8, 5),
) -> plt.Figure:
    """Histogram of GC content across all candidate guides.

    Args:
        candidates: DataFrame with 'gc_content' column.
        figsize: Figure dimensions.

    Returns:
        matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=figsize)

    if candidates.empty:
        ax.text(0.5, 0.5, "No candidates to display", ha="center", va="center")
        return fig

    ax.hist(
        candidates["gc_content"], bins=20, range=(0.3, 0.7),
        color="#26A69A", edgecolor="white", alpha=0.85,
    )

    ax.axvline(x=0.40, color="#FF5722", linestyle="--", alpha=0.7, label="Filter bounds (40–60%)")
    ax.axvline(x=0.60, color="#FF5722", linestyle="--", alpha=0.7)
    ax.axvline(x=0.50, color="#1565C0", linestyle="-", alpha=0.7, label="Ideal (50%)")

    ax.set_xlabel("GC Content")
    ax.set_ylabel("Number of Candidates")
    ax.set_title("GC Content Distribution of Guide RNA Candidates", fontweight="bold")
    ax.legend()

    plt.tight_layout()
    return fig


def plot_conservation_vs_offtarget(
    scored_df: pd.DataFrame,
    figsize: tuple = (8, 6),
) -> plt.Figure:
    """Scatter plot of conservation vs off-target hits, colored by rank.

    Args:
        scored_df: DataFrame with conservation, off_target_hits, and composite_score.
        figsize: Figure dimensions.

    Returns:
        matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=figsize)

    df = scored_df.copy()
    if df.empty:
        ax.text(0.5, 0.5, "No candidates to display", ha="center", va="center")
        return fig

    # Filter out BLAST failures
    df = df[df["off_target_hits"] >= 0]

    if df.empty:
        ax.text(0.5, 0.5, "No valid off-target data", ha="center", va="center")
        return fig

    scatter = ax.scatter(
        df["conservation"],
        df["off_target_hits"],
        c=df["composite_score"],
        cmap="RdYlGn",
        s=80,
        edgecolors="white",
        linewidth=0.5,
        alpha=0.85,
    )

    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("Composite Score")

    # Annotate top 3
    top3 = df.nlargest(3, "composite_score")
    for _, row in top3.iterrows():
        ax.annotate(
            f"#{int(row.get('rank', 0))}",
            (row["conservation"], row["off_target_hits"]),
            textcoords="offset points",
            xytext=(8, 8),
            fontsize=9,
            fontweight="bold",
            color="#1565C0",
        )

    ax.set_xlabel("Conservation Score")
    ax.set_ylabel("Off-Target Hits (Human Transcriptome)")
    ax.set_title("Conservation vs. Off-Target Specificity", fontweight="bold")

    # Ideal region annotation
    ax.annotate(
        "Ideal:\nHigh conservation\nZero off-targets",
        xy=(0.97, 0), xycoords="data",
        fontsize=8, fontstyle="italic", color="gray",
        ha="right", va="bottom",
    )

    plt.tight_layout()
    return fig
