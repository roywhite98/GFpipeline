"""Gene sequence index builder for genome-db stage.

Builds FASTA indices, gene2transcript and transcript2location mappings.
"""

from __future__ import annotations

import shutil
from collections import defaultdict
from pathlib import Path
from typing import Optional

from gfpipeline.config.schema import PipelineConfig
from gfpipeline.core.exceptions import PipelineError, StageInputError
from gfpipeline.core.file_manager import FileManager
from gfpipeline.core.logger import get_logger
from gfpipeline.core.runner import ToolRunner

log = get_logger(__name__)


def _parse_attributes(attr_str: str) -> dict[str, str]:
    """Parse GFF3 attribute column into a dict."""
    attrs: dict[str, str] = {}
    for part in attr_str.strip().split(";"):
        part = part.strip()
        if "=" in part:
            k, _, v = part.partition("=")
            attrs[k.strip()] = v.strip()
    return attrs


class GeneIndexBuilder:
    """Build gene/transcript sequence indices from GFF3 and FASTA files.

    Args:
        config: Pipeline configuration.
        runner: External tool runner (for samtools faidx).
        fm: File manager.
    """

    def __init__(self, config: PipelineConfig, runner: ToolRunner, fm: FileManager) -> None:
        self.config = config
        self.runner = runner
        self.fm = fm

    # ------------------------------------------------------------------
    # Public methods
    # ------------------------------------------------------------------

    def build_fasta_index(self, fasta_path: Path) -> None:
        """Build .fai index for a FASTA file.

        Tries samtools faidx first; falls back to BioPython SeqIO.index.

        Args:
            fasta_path: Path to the FASTA file to index.
        """
        fai_path = Path(str(fasta_path) + ".fai")
        if fai_path.exists():
            log.info("FASTA index already exists: %s", fai_path)
            return

        if shutil.which(self.config.tools.samtools):
            log.info("Building FASTA index with samtools: %s", fasta_path)
            self.runner.run([self.config.tools.samtools, "faidx", str(fasta_path)])
        else:
            log.info("samtools not found; building FASTA index with BioPython: %s", fasta_path)
            self._biopython_index(fasta_path)

    def build_gene2transcript(self, gff3_path: Path) -> dict[str, list[str]]:
        """Parse GFF3 and build gene_id -> [transcript_id, ...] mapping.

        Also writes gene2transcript.tsv to index_dir.

        Args:
            gff3_path: Path to GFF3 file.

        Returns:
            Mapping of gene_id to list of transcript IDs.
        """
        gene2transcripts: dict[str, list[str]] = defaultdict(list)

        with open(gff3_path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split("\t")
                if len(parts) < 9:
                    continue
                feat_type = parts[2]
                if feat_type not in ("mRNA", "transcript"):
                    continue
                attrs = _parse_attributes(parts[8])
                t_id = attrs.get("ID", "")
                parent = attrs.get("Parent", "")
                if t_id and parent:
                    gene2transcripts[parent].append(t_id)

        # Write TSV (is_representative will be filled by rep_index later)
        index_dir = Path(self.config.genome_db.index_dir)
        index_dir.mkdir(parents=True, exist_ok=True)
        out_path = index_dir / "gene2transcript.tsv"

        with open(out_path, "w") as fh:
            fh.write("gene_id\ttranscript_id\tis_representative\n")
            for gene_id, transcripts in sorted(gene2transcripts.items()):
                for t_id in sorted(transcripts):
                    fh.write(f"{gene_id}\t{t_id}\tFalse\n")

        log.info("Wrote gene2transcript.tsv: %s (%d genes)", out_path, len(gene2transcripts))
        return dict(gene2transcripts)

    def build_transcript2location(self, gff3_path: Path) -> None:
        """Build transcript_id -> location mapping and write transcript2location.tsv.

        Args:
            gff3_path: Path to GFF3 file.
        """
        # transcript_id -> (chrom, start, end, strand, exon_count, cds_length)
        transcripts: dict[str, dict] = {}
        exon_counts: dict[str, int] = defaultdict(int)
        cds_lengths: dict[str, int] = defaultdict(int)
        transcript_to_gene: dict[str, str] = {}

        with open(gff3_path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split("\t")
                if len(parts) < 9:
                    continue
                chrom, _, feat_type, start_s, end_s, _, strand, _, attr_s = parts[:9]
                try:
                    start = int(start_s)
                    end = int(end_s)
                except ValueError:
                    continue
                attrs = _parse_attributes(attr_s)

                if feat_type in ("mRNA", "transcript"):
                    t_id = attrs.get("ID", "")
                    parent = attrs.get("Parent", "")
                    if t_id:
                        transcripts[t_id] = {
                            "chrom": chrom,
                            "start": start,
                            "end": end,
                            "strand": strand,
                        }
                        if parent:
                            transcript_to_gene[t_id] = parent

                elif feat_type == "exon":
                    parent = attrs.get("Parent", "")
                    if parent:
                        exon_counts[parent] += 1

                elif feat_type == "CDS":
                    parent = attrs.get("Parent", "")
                    if parent:
                        cds_lengths[parent] += end - start + 1

        index_dir = Path(self.config.genome_db.index_dir)
        index_dir.mkdir(parents=True, exist_ok=True)
        out_path = index_dir / "transcript2location.tsv"

        with open(out_path, "w") as fh:
            fh.write("transcript_id\tchromosome\tstart\tend\tstrand\texon_count\tcds_length\n")
            for t_id, loc in sorted(transcripts.items()):
                fh.write(
                    f"{t_id}\t{loc['chrom']}\t{loc['start']}\t{loc['end']}\t"
                    f"{loc['strand']}\t{exon_counts.get(t_id, 0)}\t{cds_lengths.get(t_id, 0)}\n"
                )

        log.info("Wrote transcript2location.tsv: %s (%d transcripts)", out_path, len(transcripts))

    def query(self, id_: str, seq_type: str) -> str:
        """Extract a sequence by ID using the built index.

        Args:
            id_: Gene or transcript ID to query.
            seq_type: One of 'cds', 'pep', 'genome'.

        Returns:
            FASTA string for the requested sequence.

        Raises:
            PipelineError: If the ID is not found or seq_type is invalid.
        """
        fasta_map = {
            "cds": self.fm.rep_cds,
            "pep": self.fm.rep_pep,
            "genome": Path(self.config.databases.genome),
        }

        if seq_type not in fasta_map:
            raise PipelineError(f"gene-index query: unknown seq_type '{seq_type}'. Use cds/pep/genome.")

        fasta_path = fasta_map[seq_type]

        # Try samtools faidx first
        if shutil.which(self.config.tools.samtools):
            try:
                result = self.runner.run([
                    self.config.tools.samtools, "faidx", str(fasta_path), id_
                ])
                if result.stdout:
                    return result.stdout
            except Exception:
                pass

        # Fallback: BioPython
        try:
            from Bio import SeqIO
            index = SeqIO.index(str(fasta_path), "fasta")
            if id_ not in index:
                raise PipelineError(f"gene-index query: ID '{id_}' not found in {fasta_path}")
            record = index[id_]
            return f">{record.id}\n{str(record.seq)}\n"
        except ImportError:
            raise PipelineError("BioPython is required for sequence queries without samtools.")

    def run(self, force: bool = False) -> None:
        """Build all indices.

        Args:
            force: If True, rebuild even if outputs already exist.

        Raises:
            StageInputError: If required input files are missing.
        """
        db_cfg = self.config.databases
        cfg = self.config.genome_db

        # Prefer fixed GFF3 if it exists
        fixed_gff3 = Path(self.config.gff_qc.output_dir) / "gff_qc.fixed.gff3"
        gff3_path = fixed_gff3 if fixed_gff3.exists() else Path(db_cfg.gff3)

        required = {
            "GFF3": gff3_path,
            "genome": Path(db_cfg.genome),
        }
        for label, p in required.items():
            if not p.exists():
                raise StageInputError(f"gene-index: {label} file not found: {p}")

        index_dir = Path(cfg.index_dir)
        index_dir.mkdir(parents=True, exist_ok=True)

        g2t_path = index_dir / "gene2transcript.tsv"
        t2l_path = index_dir / "transcript2location.tsv"

        if not force and g2t_path.exists() and t2l_path.exists():
            log.info("Skipping gene-index (outputs already exist). Use --force to re-run.")
            return

        # Build FASTA index for genome
        self.build_fasta_index(Path(db_cfg.genome))

        # Build mapping tables
        self.build_gene2transcript(gff3_path)
        self.build_transcript2location(gff3_path)

        log.info("gene-index complete. Outputs in %s", index_dir)

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _biopython_index(self, fasta_path: Path) -> None:
        """Build a BioPython SeqIO index (creates .fai-like dict in memory)."""
        try:
            from Bio import SeqIO
            SeqIO.index(str(fasta_path), "fasta")
            log.info("BioPython index built for %s", fasta_path)
        except ImportError:
            log.warning("BioPython not available; skipping index for %s", fasta_path)
        except Exception as exc:
            log.warning("Failed to build BioPython index for %s: %s", fasta_path, exc)
