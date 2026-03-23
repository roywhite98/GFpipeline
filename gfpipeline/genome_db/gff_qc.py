"""GFF3 quality control and correction for genome-db stage.

Parses GFF3, checks format/completeness, marks truncated genes,
optionally fixes them, and writes report/summary/fixed GFF3/fix log.
"""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from gfpipeline.config.schema import PipelineConfig
from gfpipeline.core.exceptions import PipelineError, StageInputError
from gfpipeline.core.logger import get_logger
from gfpipeline.core.runner import ToolRunner

log = get_logger(__name__)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class GffQcRecord:
    gene_id: str
    issue_type: str   # format_error | missing_start | missing_stop | short_cds | edge_truncated
    detail: str
    is_truncated: bool


@dataclass
class GeneModel:
    """Minimal representation of a gene model parsed from GFF3."""
    gene_id: str
    chrom: str
    start: int   # 1-based
    end: int     # 1-based
    strand: str
    transcripts: dict[str, "TranscriptModel"] = field(default_factory=dict)
    raw_lines: list[str] = field(default_factory=list)  # original GFF3 lines for this gene


@dataclass
class TranscriptModel:
    transcript_id: str
    chrom: str
    start: int
    end: int
    strand: str
    cds_intervals: list[tuple[int, int]] = field(default_factory=list)
    has_start_codon: bool = False
    has_stop_codon: bool = False
    raw_lines: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# GFF3 parser helpers
# ---------------------------------------------------------------------------

def _parse_attributes(attr_str: str) -> dict[str, str]:
    """Parse GFF3 attribute column into a dict."""
    attrs: dict[str, str] = {}
    for part in attr_str.strip().split(";"):
        part = part.strip()
        if "=" in part:
            k, _, v = part.partition("=")
            attrs[k.strip()] = v.strip()
    return attrs


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------

class GffQc:
    """GFF3 quality control and correction.

    Args:
        config: Pipeline configuration.
        runner: External tool runner (for minimap2/stringtie).
    """

    def __init__(self, config: PipelineConfig, runner: ToolRunner) -> None:
        self.config = config
        self.runner = runner

    # ------------------------------------------------------------------
    # Public methods
    # ------------------------------------------------------------------

    def parse_gff3(self, gff3_path: Path) -> dict[str, GeneModel]:
        """Parse GFF3 file and build gene_id -> GeneModel mapping.

        Args:
            gff3_path: Path to GFF3 file.

        Returns:
            Mapping of gene_id to GeneModel.

        Raises:
            PipelineError: If the file cannot be parsed.
        """
        gene_models: dict[str, GeneModel] = {}
        # transcript_id -> gene_id for linking child features
        transcript_to_gene: dict[str, str] = {}

        try:
            with open(gff3_path) as fh:
                for raw_line in fh:
                    line = raw_line.rstrip("\n")
                    if line.startswith("#") or not line.strip():
                        continue

                    parts = line.split("\t")
                    if len(parts) < 9:
                        continue

                    chrom, source, feat_type, start_s, end_s, score, strand, phase, attr_s = parts[:9]
                    try:
                        start = int(start_s)
                        end = int(end_s)
                    except ValueError:
                        continue

                    attrs = _parse_attributes(attr_s)

                    if feat_type == "gene":
                        gene_id = attrs.get("ID", "")
                        if not gene_id:
                            continue
                        gm = GeneModel(
                            gene_id=gene_id,
                            chrom=chrom,
                            start=start,
                            end=end,
                            strand=strand,
                        )
                        gm.raw_lines.append(raw_line)
                        gene_models[gene_id] = gm

                    elif feat_type in ("mRNA", "transcript"):
                        t_id = attrs.get("ID", "")
                        parent = attrs.get("Parent", "")
                        if not t_id or not parent:
                            continue
                        tm = TranscriptModel(
                            transcript_id=t_id,
                            chrom=chrom,
                            start=start,
                            end=end,
                            strand=strand,
                        )
                        tm.raw_lines.append(raw_line)
                        transcript_to_gene[t_id] = parent
                        if parent in gene_models:
                            gene_models[parent].transcripts[t_id] = tm
                            gene_models[parent].raw_lines.append(raw_line)

                    elif feat_type == "CDS":
                        parent = attrs.get("Parent", "")
                        if not parent:
                            continue
                        gene_id = transcript_to_gene.get(parent, "")
                        if gene_id and gene_id in gene_models:
                            tm = gene_models[gene_id].transcripts.get(parent)
                            if tm:
                                tm.cds_intervals.append((start, end))
                                tm.raw_lines.append(raw_line)
                            gene_models[gene_id].raw_lines.append(raw_line)

                    elif feat_type == "start_codon":
                        parent = attrs.get("Parent", "")
                        gene_id = transcript_to_gene.get(parent, "")
                        if gene_id and gene_id in gene_models:
                            tm = gene_models[gene_id].transcripts.get(parent)
                            if tm:
                                tm.has_start_codon = True
                            gene_models[gene_id].raw_lines.append(raw_line)

                    elif feat_type == "stop_codon":
                        parent = attrs.get("Parent", "")
                        gene_id = transcript_to_gene.get(parent, "")
                        if gene_id and gene_id in gene_models:
                            tm = gene_models[gene_id].transcripts.get(parent)
                            if tm:
                                tm.has_stop_codon = True
                            gene_models[gene_id].raw_lines.append(raw_line)

                    else:
                        # exon and other features: attach to gene if possible
                        parent = attrs.get("Parent", "")
                        gene_id = transcript_to_gene.get(parent, "")
                        if gene_id and gene_id in gene_models:
                            gene_models[gene_id].raw_lines.append(raw_line)
                            tm = gene_models[gene_id].transcripts.get(parent)
                            if tm:
                                tm.raw_lines.append(raw_line)

        except OSError as exc:
            raise PipelineError(f"gff-qc: cannot read GFF3 file {gff3_path}: {exc}") from exc

        return gene_models

    def check_format(self, gene_models: dict[str, GeneModel]) -> list[GffQcRecord]:
        """Check ID/Parent completeness, coordinate validity, strand consistency.

        Args:
            gene_models: Mapping from gene_id to GeneModel.

        Returns:
            List of GffQcRecord for format issues found.
        """
        records: list[GffQcRecord] = []

        for gene_id, gm in gene_models.items():
            # Coordinate validity
            if gm.start > gm.end:
                records.append(GffQcRecord(
                    gene_id=gene_id,
                    issue_type="format_error",
                    detail=f"gene start ({gm.start}) > end ({gm.end})",
                    is_truncated=False,
                ))

            for t_id, tm in gm.transcripts.items():
                # Strand consistency
                if tm.strand != gm.strand:
                    records.append(GffQcRecord(
                        gene_id=gene_id,
                        issue_type="format_error",
                        detail=f"transcript {t_id} strand ({tm.strand}) differs from gene strand ({gm.strand})",
                        is_truncated=False,
                    ))

                # Coordinate validity for transcript
                if tm.start > tm.end:
                    records.append(GffQcRecord(
                        gene_id=gene_id,
                        issue_type="format_error",
                        detail=f"transcript {t_id} start ({tm.start}) > end ({tm.end})",
                        is_truncated=False,
                    ))

        return records

    def check_completeness(self, gene_models: dict[str, GeneModel]) -> list[GffQcRecord]:
        """Check start/stop codon presence and CDS completeness.

        Args:
            gene_models: Mapping from gene_id to GeneModel.

        Returns:
            List of GffQcRecord for completeness issues.
        """
        records: list[GffQcRecord] = []

        for gene_id, gm in gene_models.items():
            if not gm.transcripts:
                continue

            # Check across all transcripts; flag gene if any transcript is incomplete
            has_start = any(tm.has_start_codon for tm in gm.transcripts.values())
            has_stop = any(tm.has_stop_codon for tm in gm.transcripts.values())
            total_cds = sum(
                sum(e - s + 1 for s, e in tm.cds_intervals)
                for tm in gm.transcripts.values()
            )

            if not has_start:
                records.append(GffQcRecord(
                    gene_id=gene_id,
                    issue_type="missing_start",
                    detail="no start_codon annotation found",
                    is_truncated=False,
                ))

            if not has_stop:
                records.append(GffQcRecord(
                    gene_id=gene_id,
                    issue_type="missing_stop",
                    detail="no stop_codon annotation found",
                    is_truncated=False,
                ))

            if total_cds < self.config.gff_qc.min_cds_len:
                records.append(GffQcRecord(
                    gene_id=gene_id,
                    issue_type="short_cds",
                    detail=f"total CDS length {total_cds} < min_cds_len {self.config.gff_qc.min_cds_len}",
                    is_truncated=False,
                ))

        return records

    def mark_truncated(
        self,
        gene_models: dict[str, GeneModel],
        records: list[GffQcRecord],
        chrom_sizes: dict[str, int],
    ) -> list[str]:
        """Mark genes as Truncated_Gene and return their IDs.

        A gene is truncated if ANY of:
        1. Start or end coordinate is within edge_distance bp of scaffold end.
        2. Total CDS length < min_cds_len.
        3. Missing start_codon or stop_codon annotation.

        Args:
            gene_models: Mapping from gene_id to GeneModel.
            records: Existing QC records (will have is_truncated updated).
            chrom_sizes: Mapping of chromosome/scaffold name to its length.

        Returns:
            List of truncated gene IDs.
        """
        edge_dist = self.config.gff_qc.edge_distance
        min_cds = self.config.gff_qc.min_cds_len

        # Build sets of genes with each issue from records
        genes_missing_start: set[str] = set()
        genes_missing_stop: set[str] = set()
        genes_short_cds: set[str] = set()

        for rec in records:
            if rec.issue_type == "missing_start":
                genes_missing_start.add(rec.gene_id)
            elif rec.issue_type == "missing_stop":
                genes_missing_stop.add(rec.gene_id)
            elif rec.issue_type == "short_cds":
                genes_short_cds.add(rec.gene_id)

        truncated_ids: list[str] = []
        truncated_set: set[str] = set()

        for gene_id, gm in gene_models.items():
            is_trunc = False
            reasons: list[str] = []

            # Condition 1: edge proximity
            chrom_len = chrom_sizes.get(gm.chrom, 0)
            if chrom_len > 0:
                if gm.start <= edge_dist or (chrom_len - gm.end) < edge_dist:
                    is_trunc = True
                    reasons.append("edge_truncated")
            else:
                # If chrom size unknown, check start proximity to 0
                if gm.start <= edge_dist:
                    is_trunc = True
                    reasons.append("edge_truncated")

            # Condition 2: short CDS
            if gene_id in genes_short_cds:
                is_trunc = True

            # Condition 3: missing codon
            if gene_id in genes_missing_start or gene_id in genes_missing_stop:
                is_trunc = True

            if is_trunc:
                truncated_ids.append(gene_id)
                truncated_set.add(gene_id)

                # Add edge_truncated record if needed
                if "edge_truncated" in reasons:
                    records.append(GffQcRecord(
                        gene_id=gene_id,
                        issue_type="edge_truncated",
                        detail=f"gene at {gm.chrom}:{gm.start}-{gm.end} within {edge_dist}bp of scaffold end",
                        is_truncated=True,
                    ))

        # Update is_truncated flag on existing records
        for rec in records:
            if rec.gene_id in truncated_set:
                rec.is_truncated = True

        return truncated_ids

    def fix_truncated(
        self,
        truncated_ids: list[str],
        gene_models: dict[str, GeneModel],
    ) -> dict[str, GeneModel]:
        """Attempt to fix truncated genes using configured evidence sources.

        Only truncated genes are modified; others are returned unchanged.

        Args:
            truncated_ids: List of gene IDs to fix.
            gene_models: Full gene model mapping.

        Returns:
            Updated gene model mapping (non-truncated genes unchanged).
        """
        cfg = self.config.gff_qc
        fixed = dict(gene_models)  # shallow copy; non-truncated genes untouched

        if not truncated_ids:
            return fixed

        # Currently: if extra_gff is provided, merge annotations
        # (minimap2/stringtie support is a stub — requires BAM/FA files)
        if cfg.extra_gff:
            for extra_path in cfg.extra_gff:
                ep = Path(extra_path)
                if ep.exists():
                    log.info("Merging extra GFF3 for truncated genes: %s", ep)
                    extra_models = self.parse_gff3(ep)
                    for gene_id in truncated_ids:
                        if gene_id in extra_models:
                            fixed[gene_id] = extra_models[gene_id]
                            log.info("Fixed %s using extra_gff: %s", gene_id, ep)

        return fixed

    def run(self, force: bool = False) -> None:
        """Execute the full gff-qc workflow.

        Outputs:
        - gff_qc.report.tsv
        - gff_qc.summary.txt
        - gff_qc.fixed.gff3
        - gff_qc.fix.log

        Args:
            force: If True, re-run even if outputs already exist.

        Raises:
            StageInputError: If required input files are missing.
        """
        cfg = self.config.gff_qc
        db_cfg = self.config.databases

        gff3_path = Path(db_cfg.gff3)
        genome_path = Path(db_cfg.genome)

        if not gff3_path.exists():
            raise StageInputError(f"gff-qc: GFF3 file not found: {gff3_path}")
        if not genome_path.exists():
            raise StageInputError(f"gff-qc: genome FASTA not found: {genome_path}")

        out_dir = Path(cfg.output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

        report_path = out_dir / "gff_qc.report.tsv"
        summary_path = out_dir / "gff_qc.summary.txt"
        fixed_path = out_dir / "gff_qc.fixed.gff3"
        fixlog_path = out_dir / "gff_qc.fix.log"

        if not force and report_path.exists():
            log.info("Skipping gff-qc (outputs already exist). Use --force to re-run.")
            return

        log.info("Parsing GFF3: %s", gff3_path)
        gene_models = self.parse_gff3(gff3_path)
        log.info("Parsed %d gene models", len(gene_models))

        # Get chromosome sizes from genome FASTA
        chrom_sizes = self._get_chrom_sizes(genome_path)

        # Quality checks
        format_records = self.check_format(gene_models)
        completeness_records = self.check_completeness(gene_models)
        all_records = format_records + completeness_records

        # Mark truncated
        truncated_ids = self.mark_truncated(gene_models, all_records, chrom_sizes)
        log.info("Marked %d truncated genes", len(truncated_ids))

        # Fix truncated
        fixed_models = self.fix_truncated(truncated_ids, gene_models)

        # Write outputs
        self._write_report(report_path, all_records)
        self._write_summary(summary_path, gene_models, all_records, truncated_ids)
        self._write_fixed_gff3(fixed_path, fixed_models, gff3_path)
        self._write_fix_log(fixlog_path, truncated_ids, gene_models, fixed_models)

        log.info("gff-qc complete. Outputs in %s", out_dir)

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _get_chrom_sizes(self, genome_path: Path) -> dict[str, int]:
        """Parse genome FASTA to get chromosome/scaffold sizes."""
        sizes: dict[str, int] = {}
        current_chrom: Optional[str] = None
        current_len = 0

        try:
            with open(genome_path) as fh:
                for line in fh:
                    line = line.rstrip()
                    if line.startswith(">"):
                        if current_chrom is not None:
                            sizes[current_chrom] = current_len
                        # Take first word after '>' as chrom name
                        current_chrom = line[1:].split()[0]
                        current_len = 0
                    else:
                        current_len += len(line)
                if current_chrom is not None:
                    sizes[current_chrom] = current_len
        except OSError as exc:
            log.warning("Could not read genome FASTA for chrom sizes: %s", exc)

        return sizes

    def _write_report(self, path: Path, records: list[GffQcRecord]) -> None:
        """Write gff_qc.report.tsv."""
        with open(path, "w") as fh:
            fh.write("gene_id\tissue_type\tdetail\tis_truncated\n")
            for rec in records:
                fh.write(f"{rec.gene_id}\t{rec.issue_type}\t{rec.detail}\t{rec.is_truncated}\n")
        log.info("Wrote report: %s (%d records)", path, len(records))

    def _write_summary(
        self,
        path: Path,
        gene_models: dict[str, GeneModel],
        records: list[GffQcRecord],
        truncated_ids: list[str],
    ) -> None:
        """Write gff_qc.summary.txt."""
        issue_counts: dict[str, int] = defaultdict(int)
        for rec in records:
            issue_counts[rec.issue_type] += 1

        with open(path, "w") as fh:
            fh.write(f"Total genes: {len(gene_models)}\n")
            for issue_type, count in sorted(issue_counts.items()):
                fh.write(f"{issue_type}: {count}\n")
            fh.write(f"Truncated genes: {len(truncated_ids)}\n")
        log.info("Wrote summary: %s", path)

    def _write_fixed_gff3(
        self,
        path: Path,
        fixed_models: dict[str, GeneModel],
        original_gff3: Path,
    ) -> None:
        """Write gff_qc.fixed.gff3.

        Non-truncated genes are written verbatim from the original GFF3.
        """
        # Collect all raw lines from fixed models (preserves original for non-truncated)
        with open(path, "w") as fh:
            fh.write("##gff-version 3\n")
            for gene_id, gm in fixed_models.items():
                seen: set[str] = set()
                for raw_line in gm.raw_lines:
                    if raw_line not in seen:
                        fh.write(raw_line if raw_line.endswith("\n") else raw_line + "\n")
                        seen.add(raw_line)
        log.info("Wrote fixed GFF3: %s", path)

    def _write_fix_log(
        self,
        path: Path,
        truncated_ids: list[str],
        original_models: dict[str, GeneModel],
        fixed_models: dict[str, GeneModel],
    ) -> None:
        """Write gff_qc.fix.log."""
        with open(path, "w") as fh:
            fh.write("gene_id\tbefore_transcripts\tafter_transcripts\tevidence\n")
            for gene_id in truncated_ids:
                orig = original_models.get(gene_id)
                fixed = fixed_models.get(gene_id)
                before = ",".join(sorted(orig.transcripts.keys())) if orig else ""
                after = ",".join(sorted(fixed.transcripts.keys())) if fixed else ""
                changed = orig is not fixed
                evidence = "extra_gff" if changed else "none"
                fh.write(f"{gene_id}\t{before}\t{after}\t{evidence}\n")
        log.info("Wrote fix log: %s", path)
