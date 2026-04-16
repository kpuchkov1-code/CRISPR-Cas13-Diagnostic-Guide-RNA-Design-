"""Multiple sequence alignment and conservation analysis."""

import subprocess
import tempfile
import shutil
import numpy as np
from pathlib import Path
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO


def _find_muscle() -> str:
    """Locate the MUSCLE binary."""
    # Check PATH first
    muscle_path = shutil.which("muscle")
    if muscle_path:
        return muscle_path

    # Check project directory (bundled binary)
    project_root = Path(__file__).resolve().parent.parent
    for name in ("muscle.exe", "muscle"):
        candidate = project_root / name
        if candidate.exists():
            return str(candidate)

    raise FileNotFoundError(
        "MUSCLE not found. Install it and add to PATH, or place "
        "the binary in the project root directory."
    )


def run_alignment(records: list) -> object:
    """Run MUSCLE multiple sequence alignment.

    Args:
        records: List of Bio.SeqRecord objects.

    Returns:
        Bio.Align.MultipleSeqAlignment object.
    """
    muscle_bin = _find_muscle()

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False
    ) as infile:
        SeqIO.write(records, infile, "fasta")
        infile_path = infile.name

    outfile_path = infile_path.replace(".fasta", "_aligned.fasta")

    try:
        # MUSCLE v5 syntax
        result = subprocess.run(
            [muscle_bin, "-align", infile_path, "-output", outfile_path],
            capture_output=True,
            text=True,
            timeout=300,
        )
        if result.returncode != 0:
            # Fallback to MUSCLE v3 syntax
            result = subprocess.run(
                [muscle_bin, "-in", infile_path, "-out", outfile_path],
                capture_output=True,
                text=True,
                timeout=300,
            )
            if result.returncode != 0:
                raise RuntimeError(
                    f"MUSCLE alignment failed:\n{result.stderr}"
                )

        alignment = AlignIO.read(outfile_path, "fasta")
    finally:
        Path(infile_path).unlink(missing_ok=True)
        Path(outfile_path).unlink(missing_ok=True)

    print(f"Alignment complete: {len(alignment)} sequences, "
          f"{alignment.get_alignment_length()} positions")

    return alignment


def compute_conservation(alignment, window_size: int = 28) -> tuple:
    """Compute per-position conservation score using a sliding window.

    Conservation at each position is the fraction of non-gap sequences
    sharing the most common nucleotide. Positions where most sequences
    are gaps get a score of 0 (gap fraction threshold: 50%).

    Args:
        alignment: Bio.Align.MultipleSeqAlignment object.
        window_size: Size of the sliding window (matches guide RNA length).

    Returns:
        Tuple of (window_scores, position_scores):
        - window_scores: sliding window averages (length = aln_length - window_size + 1)
        - position_scores: per-position conservation (length = aln_length)
    """
    from collections import Counter

    aln_length = alignment.get_alignment_length()
    num_seqs = len(alignment)

    # Per-position conservation
    position_scores = np.zeros(aln_length)
    for i in range(aln_length):
        column = alignment[:, i].upper()
        bases = [b for b in column if b in "ACGTU"]

        # If >50% of sequences are gaps at this position, score = 0
        if len(bases) < num_seqs * 0.5:
            position_scores[i] = 0.0
            continue

        counts = Counter(bases)
        most_common_count = counts.most_common(1)[0][1]
        # Normalize by non-gap sequences (not total)
        position_scores[i] = most_common_count / len(bases)

    # Sliding window average
    if aln_length < window_size:
        return position_scores, position_scores

    window_scores = np.convolve(
        position_scores, np.ones(window_size) / window_size, mode="valid"
    )

    return window_scores, position_scores
