"""Unit tests for genome_db/gene_index.py — GeneIndexBuilder."""

from __future__ import annotations

from pathlib import Path

import pytest

from gfpipeline.config.schema import (
    DatabasesConfig,
    GenomeDbConfig,
    GffQcConfig,
    PipelineConfig,
)
from gfpipeline.core.exceptions import StageInputError
from gfpipeline.core.file_manager import FileManager
from gfpipeline.core.runner import ToolRunner
from gfpipeline.genome_db.gene_index import GeneIndexBuilder


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------

_MINIMAL_GFF3 = """\
##gff-version 3
Chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1
Chr1\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=t1.1;Parent=gene1
Chr1\t.\texon\t1000\t1500\t.\t+\t.\tID=e1;Parent=t1.1
Chr1\t.\tCDS\t1100\t1500\t.\t+\t0\tID=cds1;Parent=t1.1
Chr1\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=t1.2;Parent=gene1
Chr1\t.\texon\t1000\t2000\t.\t+\t.\tID=e2;Parent=t1.2
Chr1\t.\tCDS\t1100\t1900\t.\t+\t0\tID=cds2;Parent=t1.2
Chr2\t.\tgene\t500\t1000\t.\t-\t.\tID=gene2
Chr2\t.\tmRNA\t500\t1000\t.\t-\t.\tID=t2.1;Parent=gene2
Chr2\t.\texon\t500\t1000\t.\t-\t.\tID=e3;Parent=t2.1
Chr2\t.\tCDS\t550\t950\t.\t-\t0\tID=cds3;Parent=t2.1
"""

_MINIMAL_GENOME = ">Chr1\n" + "A" * 3000 + "\n>Chr2\n" + "T" * 2000 + "\n"


def _make_config(tmp_path: Path) -> PipelineConfig:
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
        genome_db=GenomeDbConfig(
            index_dir=str(tmp_path / "index"),
        ),
        gff_qc=GffQcConfig(output_dir=str(tmp_path / "gff_qc")),
    )


@pytest.fixture
def runner():
    return ToolRunner(dry_run=True)


# ---------------------------------------------------------------------------
# build_gene2transcript
# ---------------------------------------------------------------------------

class TestBuildGene2Transcript:
    def test_returns_correct_mapping(self, tmp_path, runner):
        config = _make_config(tmp_path)
        fm = FileManager(config)
        builder = GeneIndexBuilder(config, runner, fm)
        mapping = builder.build_gene2transcript(Path(config.databases.gff3))
        assert "gene1" in mapping
        assert set(mapping["gene1"]) == {"t1.1", "t1.2"}
        assert mapping["gene2"] == ["t2.1"]

    def test_writes_tsv_file(self, tmp_path, runner):
        config = _make_config(tmp_path)
        fm = FileManager(config)
        builder = GeneIndexBuilder(config, runner, fm)
        builder.build_gene2transcript(Path(config.databases.gff3))
        tsv = Path(config.genome_db.index_dir) / "gene2transcript.tsv"
        assert tsv.exists()

    def test_tsv_has_header(self, tmp_path, runner):
        config = _make_config(tmp_path)
        fm = FileManager(config)
        builder = GeneIndexBuilder(config, runner, fm)
        builder.build_gene2transcript(Path(config.databases.gff3))
        tsv = Path(config.genome_db.index_dir) / "gene2transcript.tsv"
        header = tsv.read_text().splitlines()[0]
        assert "gene_id" in header
        assert "transcript_id" in header

    def test_all_transcripts_in_tsv(self, tmp_path, runner):
        """All GFF3 transcripts must appear in gene2transcript.tsv (Req 16.7)."""
        config = _make_config(tmp_path)
        fm = FileManager(config)
        builder = GeneIndexBuilder(config, runner, fm)
        builder.build_gene2transcript(Path(config.databases.gff3))
        tsv = Path(config.genome_db.index_dir) / "gene2transcript.tsv"
        content = tsv.read_text()
        for t_id in ["t1.1", "t1.2", "t2.1"]:
            assert t_id in content


# ---------------------------------------------------------------------------
# build_transcript2location
# ---------------------------------------------------------------------------

class TestBuildTranscript2Location:
    def test_writes_tsv_file(self, tmp_path, runner):
        config = _make_config(tmp_path)
        fm = FileManager(config)
        builder = GeneIndexBuilder(config, runner, fm)
        builder.build_transcript2location(Path(config.databases.gff3))
        tsv = Path(config.genome_db.index_dir) / "transcript2location.tsv"
        assert tsv.exists()

    def test_tsv_has_correct_columns(self, tmp_path, runner):
        config = _make_config(tmp_path)
        fm = FileManager(config)
        builder = GeneIndexBuilder(config, runner, fm)
        builder.build_transcript2location(Path(config.databases.gff3))
        tsv = Path(config.genome_db.index_dir) / "transcript2location.tsv"
        header = tsv.read_text().splitlines()[0]
        for col in ["transcript_id", "chromosome", "start", "end", "strand", "exon_count", "cds_length"]:
            assert col in header

    def test_correct_cds_length(self, tmp_path, runner):
        config = _make_config(tmp_path)
        fm = FileManager(config)
        builder = GeneIndexBuilder(config, runner, fm)
        builder.build_transcript2location(Path(config.databases.gff3))
        tsv = Path(config.genome_db.index_dir) / "transcript2location.tsv"
        lines = tsv.read_text().splitlines()
        # t1.2 CDS: 1100-1900 = 801 bp
        t12_line = next(l for l in lines if l.startswith("t1.2"))
        cds_len = int(t12_line.split("\t")[6])
        assert cds_len == 801

    def test_correct_exon_count(self, tmp_path, runner):
        config = _make_config(tmp_path)
        fm = FileManager(config)
        builder = GeneIndexBuilder(config, runner, fm)
        builder.build_transcript2location(Path(config.databases.gff3))
        tsv = Path(config.genome_db.index_dir) / "transcript2location.tsv"
        lines = tsv.read_text().splitlines()
        t11_line = next(l for l in lines if l.startswith("t1.1"))
        exon_count = int(t11_line.split("\t")[5])
        assert exon_count == 1


# ---------------------------------------------------------------------------
# run (integration)
# ---------------------------------------------------------------------------

class TestGeneIndexBuilderRun:
    def test_run_creates_output_files(self, tmp_path, runner):
        config = _make_config(tmp_path)
        fm = FileManager(config)
        builder = GeneIndexBuilder(config, runner, fm)
        builder.run()
        index_dir = Path(config.genome_db.index_dir)
        assert (index_dir / "gene2transcript.tsv").exists()
        assert (index_dir / "transcript2location.tsv").exists()

    def test_run_raises_when_gff3_missing(self, tmp_path, runner):
        config = _make_config(tmp_path)
        Path(config.databases.gff3).unlink()
        fm = FileManager(config)
        builder = GeneIndexBuilder(config, runner, fm)
        with pytest.raises(StageInputError):
            builder.run()

    def test_run_raises_when_genome_missing(self, tmp_path, runner):
        config = _make_config(tmp_path)
        Path(config.databases.genome).unlink()
        fm = FileManager(config)
        builder = GeneIndexBuilder(config, runner, fm)
        with pytest.raises(StageInputError):
            builder.run()

    def test_run_skips_when_outputs_exist(self, tmp_path, runner):
        config = _make_config(tmp_path)
        index_dir = Path(config.genome_db.index_dir)
        index_dir.mkdir(parents=True, exist_ok=True)
        (index_dir / "gene2transcript.tsv").write_text("existing")
        (index_dir / "transcript2location.tsv").write_text("existing")
        fm = FileManager(config)
        builder = GeneIndexBuilder(config, runner, fm)
        builder.run(force=False)
        assert (index_dir / "gene2transcript.tsv").read_text() == "existing"

    def test_run_force_overwrites(self, tmp_path, runner):
        config = _make_config(tmp_path)
        index_dir = Path(config.genome_db.index_dir)
        index_dir.mkdir(parents=True, exist_ok=True)
        (index_dir / "gene2transcript.tsv").write_text("existing")
        (index_dir / "transcript2location.tsv").write_text("existing")
        fm = FileManager(config)
        builder = GeneIndexBuilder(config, runner, fm)
        builder.run(force=True)
        content = (index_dir / "gene2transcript.tsv").read_text()
        assert content != "existing"

    def test_uses_fixed_gff3_when_available(self, tmp_path, runner):
        config = _make_config(tmp_path)
        # Create a fixed GFF3 with only gene1
        fixed_gff3_dir = Path(config.gff_qc.output_dir)
        fixed_gff3_dir.mkdir(parents=True, exist_ok=True)
        fixed_gff3 = fixed_gff3_dir / "gff_qc.fixed.gff3"
        fixed_gff3.write_text(
            "##gff-version 3\n"
            "Chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1\n"
            "Chr1\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=t1.1;Parent=gene1\n"
        )
        fm = FileManager(config)
        builder = GeneIndexBuilder(config, runner, fm)
        builder.run(force=True)
        tsv = Path(config.genome_db.index_dir) / "gene2transcript.tsv"
        content = tsv.read_text()
        # Only gene1/t1.1 should be present (from fixed GFF3)
        assert "gene1" in content
        assert "gene2" not in content
