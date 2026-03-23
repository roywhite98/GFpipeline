"""Sequence extraction and parsing utilities for gfpipeline."""

from __future__ import annotations

import re
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def extract_fasta(db_path: Path, id_list: list[str], output_path: Path) -> int:
    """Extract sequences from a FASTA file by ID list, return count extracted.

    Uses BioPython SeqIO.index for efficient random access.
    """
    index = SeqIO.index(str(db_path), "fasta")
    records: list[SeqRecord] = []
    for seq_id in id_list:
        if seq_id in index:
            records.append(index[seq_id])
    SeqIO.write(records, str(output_path), "fasta")
    index.close()
    return len(records)


def transcript_to_gene_id(transcript_id: str) -> str:
    """Convert a transcript ID to a gene ID.

    Rice transcript IDs follow the pattern OsXXtYYYYYYY-NN.
    Conversion: replace 't' with 'g' and remove the '-NN' suffix.

    Examples:
        Os01t0936800-01 -> Os01g0936800
        Os12t0123456-02 -> Os12g0123456
    """
    # Remove the isoform suffix (-NN)
    gene_id = re.sub(r"-\d+$", "", transcript_id)
    # Replace the first lowercase 't' that separates chromosome from locus
    gene_id = re.sub(r"^(Os\d+)t", r"\1g", gene_id)
    return gene_id


def parse_hmm_idlist(hmm_out: Path, evalue_threshold: float) -> list[str]:
    """Parse hmmsearch --noali output and extract IDs with e-value < threshold.

    Returns a sorted, deduplicated list of sequence IDs.
    """
    ids: list[str] = []
    in_hits = False

    with open(hmm_out) as fh:
        for line in fh:
            # Detect start of the per-sequence hits table
            if line.startswith("Scores for complete sequences"):
                in_hits = True
                continue
            # End of hits table
            if in_hits and line.startswith("Domain annotation"):
                break
            if in_hits and line.startswith("------ inclusion"):
                break
            if not in_hits:
                continue
            # Skip header/separator lines
            stripped = line.strip()
            if not stripped or stripped.startswith("#") or stripped.startswith("---"):
                continue
            # Data lines: E-value  score  bias  E-value  score  bias  exp  N  Sequence ...
            parts = stripped.split()
            if len(parts) < 9:
                continue
            try:
                evalue = float(parts[0])
            except ValueError:
                continue
            if evalue < evalue_threshold:
                ids.append(parts[8])

    return sorted(set(ids))


def parse_blast_idlist(blast_out: Path) -> list[str]:
    """Parse blastp tabular output (outfmt 6) and return deduplicated subject IDs.

    Extracts column 2 (0-indexed column 1, the subject/sseqid field).
    Returns a sorted, deduplicated list.
    """
    ids: set[str] = set()
    with open(blast_out) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split("\t")
            if len(parts) >= 2:
                ids.add(parts[1])
    return sorted(ids)
