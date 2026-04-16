"""Off-target screening against the human transcriptome via NCBI BLAST."""

import time
import pandas as pd
from Bio.Blast import NCBIWWW, NCBIXML


def _blast_single_guide(sequence: str, evalue: float = 1e-5) -> list:
    """BLAST a single guide sequence against human RefSeq RNA.

    Args:
        sequence: Guide RNA sequence (will be converted to DNA for BLAST).
        evalue: E-value threshold.

    Returns:
        List of dicts with hit details.
    """
    # Convert RNA to DNA for BLAST
    dna_seq = sequence.upper().replace("U", "T")

    result_handle = NCBIWWW.qblast(
        program="blastn",
        database="refseq_rna",
        sequence=dna_seq,
        entrez_query="Homo sapiens[orgn]",
        expect=evalue,
        word_size=7,
        megablast=False,
    )

    blast_records = NCBIXML.parse(result_handle)
    record = next(blast_records)
    result_handle.close()

    hits = []
    for alignment in record.alignments:
        for hsp in alignment.hsps:
            identity = hsp.identities / hsp.align_length
            if identity >= 0.80 and hsp.align_length >= 20:
                hits.append({
                    "target": alignment.title[:100],
                    "identity": round(identity, 3),
                    "align_length": hsp.align_length,
                    "mismatches": hsp.align_length - hsp.identities,
                    "evalue": hsp.expect,
                })

    return hits


def screen_off_targets(
    candidates: pd.DataFrame,
    evalue: float = 1e-5,
    delay: float = 3.0,
    max_guides: int = 20,
) -> pd.DataFrame:
    """Screen candidate guides for off-target hits in the human transcriptome.

    Args:
        candidates: DataFrame from extract_candidates.
        evalue: BLAST e-value cutoff.
        delay: Seconds to wait between BLAST queries (NCBI rate limit).
        max_guides: Maximum number of top candidates to screen (BLAST is slow).

    Returns:
        Updated DataFrame with off_target_hits count and off_target_details columns.
    """
    if candidates.empty:
        return candidates

    # Only screen top candidates by conservation (BLAST is slow)
    df = candidates.head(max_guides).copy()

    hit_counts = []
    hit_details = []

    total = len(df)
    for idx, (_, row) in enumerate(df.iterrows()):
        print(f"  BLASTing guide {idx + 1}/{total}: {row['sequence'][:20]}...")

        try:
            hits = _blast_single_guide(row["sequence"], evalue=evalue)
            hit_counts.append(len(hits))
            hit_details.append(hits)
        except Exception as e:
            print(f"    BLAST failed for guide at position {row['position']}: {e}")
            hit_counts.append(-1)  # -1 indicates BLAST failure
            hit_details.append([])

        if idx < total - 1:
            time.sleep(delay)

    df["off_target_hits"] = hit_counts
    df["off_target_details"] = hit_details

    passed = (df["off_target_hits"] == 0).sum()
    flagged = (df["off_target_hits"] > 0).sum()
    failed = (df["off_target_hits"] == -1).sum()

    print(f"\nOff-target screening complete:")
    print(f"  {passed} guides with no human hits")
    print(f"  {flagged} guides flagged with off-target hits")
    if failed > 0:
        print(f"  {failed} guides where BLAST query failed")

    return df
