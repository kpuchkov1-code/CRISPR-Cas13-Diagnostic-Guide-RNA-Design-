# CRISPR-Cas13 Diagnostic Guide RNA Design Tool

Automated pipeline for designing CRISPR-Cas13 guide RNAs (crRNAs) targeting conserved pathogen genomic regions for nucleic acid diagnostics.

## Pipeline

1. **Sequence Retrieval** -- Fetches pathogen sequences from NCBI Entrez
2. **Multiple Sequence Alignment** -- MUSCLE alignment for positional conservation
3. **Guide RNA Extraction** -- 28-nt candidates from conserved regions with biophysical filtering
4. **Off-Target Screening** -- BLAST against human RefSeq transcriptome
5. **Multi-Criteria Scoring** -- Weighted composite ranking (conservation, specificity, GC, homopolymer, self-folding)

## Requirements

- Python 3.9+
- [MUSCLE](https://github.com/rcedgar/muscle) (v3 or v5, must be on PATH)
- NCBI internet access (Entrez + BLAST)

## Setup

```bash
pip install -r requirements.txt
```

Verify MUSCLE is installed:

```bash
muscle -version
```

## Usage

Open the Jupyter notebook:

```bash
cd notebooks
jupyter notebook cas13_guide_design.ipynb
```

Configure the target pathogen in the notebook:

```python
PATHOGEN = "Zika virus"     # Any pathogen name or NCBI taxonomy ID
GENE_TARGET = "NS5"         # Optional: specific gene
```

Run all cells. The pipeline will retrieve sequences, align them, design guide RNAs, screen for off-targets, and output ranked candidates with visualizations.

## Output

- Ranked table of guide RNA candidates with composite scores
- Conservation heatmap across the pathogen genome
- Score breakdown chart for top guides
- GC content distribution
- Conservation vs. off-target specificity scatter plot
- CSV export of results

## Project Structure

```
cas13-guide-designer/
├── cas13_design/
│   ├── __init__.py          # Package exports
│   ├── retrieval.py         # NCBI sequence fetching
│   ├── alignment.py         # MUSCLE MSA + conservation
│   ├── guide_design.py      # Candidate extraction & filtering
│   ├── off_target.py        # BLAST off-target screening
│   ├── scoring.py           # Multi-criteria ranking
│   └── visualization.py     # Plots and figures
├── notebooks/
│   └── cas13_guide_design.ipynb
├── requirements.txt
└── README.md
```
