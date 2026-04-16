"""Guide RNA candidate extraction and filtering."""

import numpy as np
import pandas as pd
from Bio.Align import MultipleSeqAlignment
from collections import Counter


def _consensus_at(alignment, start: int, length: int) -> str:
    """Extract consensus sequence from alignment at given region."""
    consensus = []
    for i in range(start, start + length):
        column = alignment[:, i].upper()
        bases = [b for b in column if b in "ACGTU"]
        if not bases:
            consensus.append("N")
            continue
        counts = Counter(bases)
        consensus.append(counts.most_common(1)[0][0])
    return "".join(consensus)


def _gc_content(seq: str) -> float:
    """Calculate GC content of a sequence."""
    seq = seq.upper().replace("U", "T")
    gc = sum(1 for b in seq if b in "GC")
    return gc / len(seq) if seq else 0.0


def _max_homopolymer(seq: str) -> int:
    """Find the longest homopolymer run in a sequence."""
    if not seq:
        return 0
    max_run = 1
    current_run = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            current_run += 1
            max_run = max(max_run, current_run)
        else:
            current_run = 1
    return max_run


def _self_complementarity(seq: str) -> float:
    """Score self-complementarity of a guide RNA.

    Counts the number of positions where the guide could form
    internal base pairs (comparing the first half with the
    reverse complement of the second half).

    Returns a score from 0 (no self-complementarity) to 1 (fully self-complementary).
    """
    complement = {"A": "U", "U": "A", "G": "C", "C": "G", "T": "A"}
    seq = seq.upper().replace("T", "U")
    n = len(seq)
    half = n // 2

    first_half = seq[:half]
    second_half_rc = "".join(
        complement.get(b, "N") for b in reversed(seq[half:2 * half])
    )

    matches = sum(1 for a, b in zip(first_half, second_half_rc) if a == b)
    return matches / half if half > 0 else 0.0


def extract_candidates(
    alignment,
    window_scores: np.ndarray,
    guide_length: int = 28,
    conservation_threshold: float = 0.90,
    gc_min: float = 0.40,
    gc_max: float = 0.60,
    max_homopolymer: int = 4,
) -> pd.DataFrame:
    """Extract and filter guide RNA candidates from conserved regions.

    Args:
        alignment: Bio.Align.MultipleSeqAlignment object.
        window_scores: Conservation scores from compute_conservation.
        guide_length: Length of guide RNA (default 28nt for Cas13a).
        conservation_threshold: Minimum mean conservation to consider a region.
        gc_min: Minimum acceptable GC content.
        gc_max: Maximum acceptable GC content.
        max_homopolymer: Maximum allowed homopolymer run length.

    Returns:
        DataFrame with columns: sequence, position, conservation,
        gc_content, homopolymer_run, self_comp.
    """
    candidates = []

    for i, score in enumerate(window_scores):
        if score < conservation_threshold:
            continue

        seq = _consensus_at(alignment, i, guide_length)

        # Skip sequences with ambiguous bases or gaps
        if "N" in seq or "-" in seq:
            continue

        gc = _gc_content(seq)
        if gc < gc_min or gc > gc_max:
            continue

        homo = _max_homopolymer(seq)
        if homo > max_homopolymer:
            continue

        self_comp = _self_complementarity(seq)

        candidates.append({
            "sequence": seq,
            "position": i,
            "conservation": round(score, 4),
            "gc_content": round(gc, 4),
            "homopolymer_run": homo,
            "self_comp": round(self_comp, 4),
        })

    df = pd.DataFrame(candidates)

    if df.empty:
        print("No candidates passed all filters. Consider relaxing thresholds.")
        return df

    # Remove duplicate sequences (keep the one with highest conservation)
    df = (
        df.sort_values("conservation", ascending=False)
        .drop_duplicates(subset="sequence", keep="first")
        .reset_index(drop=True)
    )

    print(f"Extracted {len(df)} unique guide RNA candidates")
    return df
