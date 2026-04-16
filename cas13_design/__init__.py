"""CRISPR-Cas13 Diagnostic Guide RNA Design Tool.

Automated pipeline for designing guide RNAs targeting conserved
pathogen genomic regions for nucleic acid diagnostics.
"""

from cas13_design.retrieval import fetch_sequences
from cas13_design.alignment import run_alignment, compute_conservation
from cas13_design.guide_design import extract_candidates
from cas13_design.off_target import screen_off_targets
from cas13_design.scoring import score_guides, rank_guides
from cas13_design.visualization import (
    plot_conservation_heatmap,
    plot_score_breakdown,
    plot_gc_distribution,
    plot_conservation_vs_offtarget,
)

__version__ = "0.1.0"
