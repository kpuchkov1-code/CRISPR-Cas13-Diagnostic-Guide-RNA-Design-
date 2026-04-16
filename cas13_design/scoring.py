"""Multi-criteria guide RNA scoring and ranking."""

import numpy as np
import pandas as pd


def score_guides(
    candidates: pd.DataFrame,
    w_conservation: float = 0.30,
    w_offtarget: float = 0.25,
    w_gc: float = 0.15,
    w_homopolymer: float = 0.15,
    w_selfcomp: float = 0.15,
) -> pd.DataFrame:
    """Compute weighted composite scores for guide RNA candidates.

    Scoring criteria:
        - Conservation (30%): Higher is better.
        - Off-target safety (25%): Fewer hits is better.
        - GC optimality (15%): Closer to 50% is better.
        - Homopolymer penalty (15%): Shorter max run is better.
        - Self-complementarity penalty (15%): Lower is better.

    Args:
        candidates: DataFrame with guide data including off_target_hits.
        w_conservation: Weight for conservation score.
        w_offtarget: Weight for off-target safety.
        w_gc: Weight for GC content optimality.
        w_homopolymer: Weight for homopolymer penalty.
        w_selfcomp: Weight for self-complementarity penalty.

    Returns:
        DataFrame with individual component scores and composite_score.
    """
    df = candidates.copy()

    if df.empty:
        return df

    # 1. Conservation score (already 0-1, higher is better)
    df["score_conservation"] = df["conservation"]

    # 2. Off-target safety (invert: 0 hits → 1.0, more hits → lower score)
    max_hits = df["off_target_hits"].replace(-1, np.nan).max()
    if pd.isna(max_hits) or max_hits == 0:
        df["score_offtarget"] = df["off_target_hits"].apply(
            lambda x: 1.0 if x == 0 else (0.5 if x == -1 else 0.0)
        )
    else:
        df["score_offtarget"] = df["off_target_hits"].apply(
            lambda x: 0.5 if x == -1 else 1.0 - (x / (max_hits + 1))
        )

    # 3. GC optimality (distance from ideal 50%, mapped to 0-1)
    df["score_gc"] = 1.0 - (abs(df["gc_content"] - 0.50) / 0.50)
    df["score_gc"] = df["score_gc"].clip(0, 1)

    # 4. Homopolymer penalty (shorter runs → higher score)
    # Max possible in a 28-mer is 28; practical max for filtered guides is 4
    df["score_homopolymer"] = 1.0 - ((df["homopolymer_run"] - 1) / 3)
    df["score_homopolymer"] = df["score_homopolymer"].clip(0, 1)

    # 5. Self-complementarity penalty (lower self-comp → higher score)
    df["score_selfcomp"] = 1.0 - df["self_comp"]

    # Composite weighted score
    df["composite_score"] = (
        w_conservation * df["score_conservation"]
        + w_offtarget * df["score_offtarget"]
        + w_gc * df["score_gc"]
        + w_homopolymer * df["score_homopolymer"]
        + w_selfcomp * df["score_selfcomp"]
    )

    return df


def rank_guides(scored_df: pd.DataFrame, top_n: int = 10) -> pd.DataFrame:
    """Rank guides by composite score and return top candidates.

    Args:
        scored_df: DataFrame from score_guides.
        top_n: Number of top guides to return.

    Returns:
        Sorted DataFrame with rank column.
    """
    df = (
        scored_df.sort_values("composite_score", ascending=False)
        .head(top_n)
        .reset_index(drop=True)
    )
    df["rank"] = range(1, len(df) + 1)

    display_cols = [
        "rank", "sequence", "position", "composite_score",
        "score_conservation", "score_offtarget", "score_gc",
        "score_homopolymer", "score_selfcomp",
        "conservation", "gc_content", "off_target_hits",
    ]
    existing_cols = [c for c in display_cols if c in df.columns]

    return df[existing_cols]
