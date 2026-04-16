# CRISPR-Cas13 Diagnostic Guide RNA Design Tool — Design Spec

**Date:** 2026-04-16
**Status:** Approved

## Overview

A Python pipeline for automated design of CRISPR-Cas13 guide RNAs targeting conserved pathogen genomic regions for nucleic acid diagnostics. User specifies any pathogen; the tool retrieves sequences from NCBI, aligns them, extracts conserved regions, screens for off-target human hits, and ranks candidate guides by a multi-criteria composite score.

## Interface

Jupyter notebook — modular Python package (`cas13_design/`) handles logic, notebook serves as the narrative presentation layer.

## Pipeline Stages

1. **Retrieval** — Entrez esearch/efetch, up to 50 sequences, user provides pathogen name/taxID + optional gene target
2. **Alignment** — MUSCLE via subprocess, per-position conservation via sliding window (28nt)
3. **Guide Design** — Extract 28-mers from regions with ≥90% conservation, filter GC 40-60%, no homopolymer runs ≥5nt, self-complementarity penalty
4. **Off-Target Screening** — Remote BLAST vs human RefSeq RNA (refseq_rna, Homo sapiens), flag guides with ≥80% identity over ≥20nt
5. **Scoring & Ranking** — Weighted composite: conservation 30%, off-target safety 25%, GC optimality 15%, homopolymer penalty 15%, self-folding penalty 15%

## Package Modules

- `retrieval.py` — NCBI sequence fetching
- `alignment.py` — MUSCLE MSA + conservation scoring
- `guide_design.py` — Candidate extraction, filtering, self-complementarity
- `off_target.py` — BLAST against human transcriptome
- `scoring.py` — Multi-criteria ranking
- `visualization.py` — Conservation heatmap, score bar chart, GC distribution, conservation vs off-target scatter

## Notebook Sections

1. Introduction (what is Cas13 diagnostics)
2. Setup & Configuration
3. Sequence Retrieval
4. Multiple Sequence Alignment
5. Guide RNA Candidate Generation
6. Off-Target Screening
7. Scoring & Ranking
8. Visualization
9. Results Summary (top 5 guides, CSV export)

## Dependencies

biopython, pandas, numpy, matplotlib, seaborn, MUSCLE (external binary on PATH)

## Project Structure

```
cas13-guide-designer/
├── cas13_design/
│   ├── __init__.py
│   ├── retrieval.py
│   ├── alignment.py
│   ├── guide_design.py
│   ├── off_target.py
│   ├── scoring.py
│   └── visualization.py
├── notebooks/
│   └── cas13_guide_design.ipynb
├── requirements.txt
└── README.md
```
