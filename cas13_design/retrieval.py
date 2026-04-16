"""Sequence retrieval from NCBI via Entrez."""

import statistics
from typing import Optional
from Bio import Entrez, SeqIO
from io import StringIO


Entrez.email = "kp425@ic.ac.uk"


def _filter_by_length(records: list, tolerance: float = 0.5) -> list:
    """Keep sequences within a tolerance band around the median length.

    This removes outliers (e.g., fragments or full genomes mixed with genes)
    that would produce gappy, low-quality alignments.

    Args:
        records: List of SeqRecord objects.
        tolerance: Fraction of median length to use as the band (0.5 = +/-50%).

    Returns:
        Filtered list of SeqRecord objects.
    """
    if len(records) <= 3:
        return records

    lengths = [len(r.seq) for r in records]
    median_len = statistics.median(lengths)
    low = median_len * (1 - tolerance)
    high = median_len * (1 + tolerance)

    filtered = [r for r in records if low <= len(r.seq) <= high]

    # Fall back to original if filtering is too aggressive
    if len(filtered) < 3:
        return records

    return filtered


def fetch_sequences(
    pathogen: str,
    gene: Optional[str] = None,
    max_seqs: int = 50,
    db: str = "nucleotide",
) -> list:
    """Fetch pathogen sequences from NCBI.

    Retrieves sequences and filters for length uniformity to ensure
    high-quality alignment. Excludes patent sequences.

    Args:
        pathogen: Pathogen name or NCBI taxonomy ID.
        gene: Optional gene name to narrow the search (e.g., "N gene", "RdRp").
        max_seqs: Maximum number of sequences to retrieve.
        db: NCBI database to search.

    Returns:
        List of Bio.SeqRecord objects.
    """
    query = f"{pathogen}[Organism]"
    if gene:
        query += f" AND {gene}[Gene Name]"
    query += " AND 200:10000[Sequence Length]"
    query += " NOT patent[Title]"

    # Request more than needed so we can filter
    fetch_count = min(max_seqs * 3, 150)

    handle = Entrez.esearch(db=db, term=query, retmax=fetch_count, usehistory="y")
    search_results = Entrez.read(handle)
    handle.close()

    count = int(search_results["Count"])
    if count == 0:
        raise ValueError(
            f"No sequences found for query: {query}. "
            "Try a different pathogen name or gene target."
        )

    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]

    fetch_handle = Entrez.efetch(
        db=db,
        rettype="fasta",
        retmode="text",
        retmax=fetch_count,
        webenv=webenv,
        query_key=query_key,
    )
    fasta_data = fetch_handle.read()
    fetch_handle.close()

    records = list(SeqIO.parse(StringIO(fasta_data), "fasta"))

    # Filter for length uniformity
    records = _filter_by_length(records)

    # Cap at max_seqs
    records = records[:max_seqs]

    if len(records) < 3:
        raise ValueError(
            f"Only {len(records)} sequences retrieved — need at least 3 for "
            "meaningful alignment. Try broadening your search."
        )

    lengths = [len(r.seq) for r in records]
    print(f"Retrieved {len(records)} sequences (query: {query})")
    print(f"Length range: {min(lengths)}–{max(lengths)} bp "
          f"(median {statistics.median(lengths):.0f})")
    for rec in records[:5]:
        print(f"  {rec.id}: {len(rec.seq)} bp — {rec.description[:80]}")
    if len(records) > 5:
        print(f"  ... and {len(records) - 5} more")

    return records
