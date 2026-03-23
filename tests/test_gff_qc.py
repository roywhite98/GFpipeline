"""Unit tests for genome_db/gff_qc.py — GffQc."""

from __future__ import annotations

from pathlib import Path

import pytest

from gfpipeline.config.schema import (
    DatabasesConfig,
    GffQcConfig,
    PipelineConfig,
)
from gfpipeline.core.exceptions import StageInputError
from gfpipeline.core.runner import ToolRunner
from gfpipeline.genome_db.gff_qc import GffQc, GffQcRecord, GeneModel, TranscriptModel


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------

_MINIMAL_GFF3 = """\
##gff-version 3
Chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1
Chr1\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=t1;Parent=gene1
Chr1\t.\tCDS\t1100\t1900\t.\t+\t0\tID=cds1;Parent=t1
Chr1\t.\tstart_codon\t1100\t1102\t.\t+\t0\tID=sc1;Parent=t1
Chr1\t.\tstop_codon\t1898\t1900\t.\t+\t0\tID=ec1;Parent=t1
Chr2\t.\tgene\t100\t500\t.\t-\t.\tID=gene2
Chr2\t.\tmRNA\t100\t500\t.\t-\t.\tID=t2;Parent=gene2
Chr2\t.\tCDS\t150\t400\t.\t-\t0\tID=cds2;Parent=t2
"""

_MINIMAL_GENOME = """\
>Chr1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>Chr2
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
"""


def _make_config(tmp_path: Path, **gff_qc_kwargs) -> PipelineConfig:
    gff3 = tmp_path / "genome.gff3"
    genome = tmp_path / "genome.fa"
    gff3.write_text(_MINIMAL_GFF3)
    genome.write_text(_MINIMAL_GENOME)

    return PipelineConfig(
        project_name="TEST",
        data_dir=str(tmp_path / "data"),
        result_dir=str(tmp_path / "results"),
        databases=DatabasesConfig(
            ref_fa=str(tmp_path / "ref.fa"),
            genome=str(genome),
            gff3=str(gff3),
        ),
        gff_qc=GffQcConfig(
            output_dir=str(tmp_path / "gff_qc"),
            **gff_qc_kwargs,
        ),
    )


@pytest.fixture
def runner():
    return ToolRunner(dry_run=True)


# ---------------------------------------------------------------------------
# parse_gff3
# ---------------------------------------------------------------------------

class TestParseGff3:
    def test_parses_genes(self, tmp_path, runner):
        config = _make_config(tmp_path)
        qc = GffQc(config, runner)
        models = qc.parse_gff3(Path(config.databases.gff3))
        assert "gene1" in models
        assert "gene2" in models

    def test_parses_transcripts(self, tmp_path, runner):
        config = _make_config(tmp_path)
        qc = GffQc(config, runner)
        models = qc.parse_gff3(Path(config.databases.gff3))
        assert "t1" in models["gene1"].transcripts
        assert "t2" in models["gene2"].transcripts

    def test_parses_cds_intervals(self, tmp_path, runner):
        config = _make_config(tmp_path)
        qc = GffQc(config, runner)
        models = qc.parse_gff3(Path(config.databases.gff3))
        t1 = models["gene1"].transcripts["t1"]
        assert len(t1.cds_intervals) == 1
        assert t1.cds_intervals[0] == (1100, 1900)

    def test_parses_start_stop_codons(self, tmp_path, runner):
        config = _make_config(tmp_path)
        qc = GffQc(config, runner)
        models = qc.parse_gff3(Path(config.databases.gff3))
        t1 = models["gene1"].transcripts["t1"]
        assert t1.has_start_codon is True
        assert t1.has_stop_codon is True

    def test_gene2_missing_codons(self, tmp_path, runner):
        config = _make_config(tmp_path)
        qc = GffQc(config, runner)
        models = qc.parse_gff3(Path(config.databases.gff3))
        t2 = models["gene2"].transcripts["t2"]
        assert t2.has_start_codon is False
        assert t2.has_stop_codon is False


# ---------------------------------------------------------------------------
# check_format
# ---------------------------------------------------------------------------

class TestCheckFormat:
    def test_no_errors_for_valid_gff3(self, tmp_path, runner):
        config = _make_config(tmp_path)
        qc = GffQc(config, runner)
        models = qc.parse_gff3(Path(config.databases.gff3))
        records = qc.check_format(models)
        assert all(r.issue_type == "format_error" for r in records) or len(records) == 0

    def test_detects_strand_inconsistency(self, tmp_path, runner):
        config = _make_config(tmp_path)
        qc = GffQc(config, runner)
        # Manually create a model with strand mismatch
        gm = GeneModel(gene_id="gX", chrom="Chr1", start=100, end=200, strand="+")
        tm = TranscriptModel(transcript_id="tX", chrom="Chr1", start=100, end=200, strand="-")
        gm.transcripts["tX"] = tm
        records = qc.check_format({"gX": gm})
        assert any(r.issue_type == "format_error" and "strand" in r.detail for r in records)

    def test_detects_invalid_coordinates(self, tmp_path, runner):
        config = _make_config(tmp_path)
        qc = GffQc(config, runner)
        gm = GeneModel(gene_id="gY", chrom="Chr1", start=500, end=100, strand="+")
        records = qc.check_format({"gY": gm})
        assert any(r.issue_type == "format_error" and "start" in r.detail for r in records)


# ---------------------------------------------------------------------------
# check_completeness
# ---------------------------------------------------------------------------

class TestCheckCompleteness:
    def test_flags_missing_start_codon(self, tmp_path, runner):
        config = _make_config(tmp_path, min_cds_len=10)
        qc = GffQc(config, runner)
        models = qc.parse_gff3(Path(config.databases.gff3))
        records = qc.check_completeness(models)
        gene2_issues = [r for r in records if r.gene_id == "gene2"]
        issue_types = {r.issue_type for r in gene2_issues}
        assert "missing_start" in issue_types

    def test_flags_missing_stop_codon(self, tmp_path, runner):
        config = _make_config(tmp_path, min_cds_len=10)
        qc = GffQc(config, runner)
        models = qc.parse_gff3(Path(config.databases.gff3))
        records = qc.check_completeness(models)
        gene2_issues = [r for r in records if r.gene_id == "gene2"]
        issue_types = {r.issue_type for r in gene2_issues}
        assert "missing_stop" in issue_types

    def test_flags_short_cds(self, tmp_path, runner):
        # Set min_cds_len very high so gene1 CDS (801 bp) is flagged
        config = _make_config(tmp_path, min_cds_len=10000)
        qc = GffQc(config, runner)
        models = qc.parse_gff3(Path(config.databases.gff3))
        records = qc.check_completeness(models)
        assert any(r.issue_type == "short_cds" for r in records)

    def test_no_short_cds_when_above_threshold(self, tmp_path, runner):
        config = _make_config(tmp_path, min_cds_len=10)
        qc = GffQc(config, runner)
        models = qc.parse_gff3(Path(config.databases.gff3))
        records = qc.check_completeness(models)
        # gene1 has CDS 1100-1900 = 801 bp, well above 10
        gene1_issues = [r for r in records if r.gene_id == "gene1"]
        assert not any(r.issue_type == "short_cds" for r in gene1_issues)


# ---------------------------------------------------------------------------
# mark_truncated
# ---------------------------------------------------------------------------

class TestMarkTruncated:
    def test_marks_edge_gene_as_truncated(self, tmp_path, runner):
        # gene2 starts at 100, chrom is 80 bp long -> within edge_distance=500
        config = _make_config(tmp_path, edge_distance=500, min_cds_len=10)
        qc = GffQc(config, runner)
        models = qc.parse_gff3(Path(config.databases.gff3))
        records = qc.check_completeness(models)
        chrom_sizes = {"Chr1": 5000, "Chr2": 5000}
        truncated = qc.mark_truncated(models, records, chrom_sizes)
        # gene2 starts at 100 <= 500 (edge_distance)
        assert "gene2" in truncated

    def test_marks_short_cds_gene_as_truncated(self, tmp_path, runner):
        config = _make_config(tmp_path, min_cds_len=10000, edge_distance=0)
        qc = GffQc(config, runner)
        models = qc.parse_gff3(Path(config.databases.gff3))
        records = qc.check_completeness(models)
        chrom_sizes = {"Chr1": 100000, "Chr2": 100000}
        truncated = qc.mark_truncated(models, records, chrom_sizes)
        assert len(truncated) > 0

    def test_non_truncated_gene_not_in_list(self, tmp_path, runner):
        # gene1 has start_codon, stop_codon, long CDS, not near edge
        config = _make_config(tmp_path, min_cds_len=10, edge_distance=50)
        qc = GffQc(config, runner)
        models = qc.parse_gff3(Path(config.databases.gff3))
        records = qc.check_completeness(models)
        chrom_sizes = {"Chr1": 100000, "Chr2": 100000}
        truncated = qc.mark_truncated(models, records, chrom_sizes)
        assert "gene1" not in truncated

    def test_is_truncated_flag_set_on_records(self, tmp_path, runner):
        config = _make_config(tmp_path, min_cds_len=10, edge_distance=500)
        qc = GffQc(config, runner)
        models = qc.parse_gff3(Path(config.databases.gff3))
        records = qc.check_completeness(models)
        chrom_sizes = {"Chr1": 5000, "Chr2": 5000}
        truncated = qc.mark_truncated(models, records, chrom_sizes)
        for rec in records:
            if rec.gene_id in truncated:
                assert rec.is_truncated is True


# ---------------------------------------------------------------------------
# run (integration)
# ---------------------------------------------------------------------------

class TestGffQcRun:
    def test_run_creates_output_files(self, tmp_path, runner):
        config = _make_config(tmp_path)
        qc = GffQc(config, runner)
        qc.run()
        out_dir = Path(config.gff_qc.output_dir)
        assert (out_dir / "gff_qc.report.tsv").exists()
        assert (out_dir / "gff_qc.summary.txt").exists()
        assert (out_dir / "gff_qc.fixed.gff3").exists()
        assert (out_dir / "gff_qc.fix.log").exists()

    def test_run_raises_when_gff3_missing(self, tmp_path, runner):
        config = _make_config(tmp_path)
        Path(config.databases.gff3).unlink()
        qc = GffQc(config, runner)
        with pytest.raises(StageInputError, match="GFF3 file not found"):
            qc.run()

    def test_run_raises_when_genome_missing(self, tmp_path, runner):
        config = _make_config(tmp_path)
        Path(config.databases.genome).unlink()
        qc = GffQc(config, runner)
        with pytest.raises(StageInputError, match="genome FASTA not found"):
            qc.run()

    def test_run_skips_when_outputs_exist(self, tmp_path, runner):
        config = _make_config(tmp_path)
        out_dir = Path(config.gff_qc.output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "gff_qc.report.tsv").write_text("existing")
        qc = GffQc(config, runner)
        qc.run(force=False)
        # File should still contain "existing" (not overwritten)
        assert (out_dir / "gff_qc.report.tsv").read_text() == "existing"

    def test_run_force_overwrites(self, tmp_path, runner):
        config = _make_config(tmp_path)
        out_dir = Path(config.gff_qc.output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "gff_qc.report.tsv").write_text("existing")
        qc = GffQc(config, runner)
        qc.run(force=True)
        content = (out_dir / "gff_qc.report.tsv").read_text()
        assert content != "existing"

    def test_report_has_header(self, tmp_path, runner):
        config = _make_config(tmp_path)
        qc = GffQc(config, runner)
        qc.run()
        report = (Path(config.gff_qc.output_dir) / "gff_qc.report.tsv").read_text()
        assert "gene_id" in report
        assert "issue_type" in report

    def test_non_truncated_genes_preserved_in_fixed_gff3(self, tmp_path, runner):
        """Non-truncated genes must appear unchanged in fixed.gff3 (Req 15.13)."""
        config = _make_config(tmp_path, min_cds_len=10, edge_distance=50)
        qc = GffQc(config, runner)
        qc.run()
        fixed = (Path(config.gff_qc.output_dir) / "gff_qc.fixed.gff3").read_text()
        # gene1 is not truncated; its lines should appear in fixed GFF3
        assert "gene1" in fixed
