"""Unit tests for genome_db/rep_index.py — RepIndexBuilder."""

from __future__ import annotations

from pathlib import Path

import pytest

from gfpipeline.config.schema import (
    DatabasesConfig,
    GenomeDbConfig,
    PipelineConfig,
)
from gfpipeline.core.exceptions import StageInputError
from gfpipeline.core.file_manager import FileManager
from gfpipeline.genome_db.rep_index import RepIndexBuilder


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------

# GFF3: gene1 has two transcripts (t1.1 shorter CDS, t1.2 longer CDS)
#        gene2 has one transcript on minus strand
#        gene3 has one transcript but NO CDS (should be skipped)
_MINIMAL_GFF3 = (
    "##gff-version 3\n"
    "Chr1\t.\tgene\t1\t30\t.\t+\t.\tID=gene1\n"
    "Chr1\t.\tmRNA\t1\t30\t.\t+\t.\tID=t1.1;Parent=gene1\n"
    "Chr1\t.\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=t1.1\n"
    "Chr1\t.\tmRNA\t1\t30\t.\t+\t.\tID=t1.2;Parent=gene1\n"
    "Chr1\t.\tCDS\t1\t15\t.\t+\t0\tID=cds2;Parent=t1.2\n"
    "Chr2\t.\tgene\t1\t9\t.\t-\t.\tID=gene2\n"
    "Chr2\t.\tmRNA\t1\t9\t.\t-\t.\tID=t2.1;Parent=gene2\n"
    "Chr2\t.\tCDS\t1\t9\t.\t-\t0\tID=cds3;Parent=t2.1\n"
    "Chr3\t.\tgene\t1\t20\t.\t+\t.\tID=gene3\n"
    "Chr3\t.\tmRNA\t1\t20\t.\t+\t.\tID=t3.1;Parent=gene3\n"
)

# Genome sequences — must contain valid codons for translation tests
_MINIMAL_GENOME = (
    ">Chr1\n"
    "ATGCCCGGGAAATTTCCCGGGAAATTTCCC\n"
    ">Chr2\n"
    "ATGAAACCC\n"
    ">Chr3\n"
    "ATGAAATTTGGGCCCAAATTT\n"
)


def _make_config(
    tmp_path: Path,
    rep_selection: str = "longest_cds",
    genome_name: str | None = None,
) -> PipelineConfig:
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
            rep_index_dir=str(tmp_path / "rep_index"),
            rep_selection=rep_selection,
            genome_name=genome_name,
        ),
    )


# ---------------------------------------------------------------------------
# select_representative
# ---------------------------------------------------------------------------

class TestSelectRepresentative:
    def setup_method(self):
        self._config = PipelineConfig(
            project_name="X",
            data_dir="data",
            result_dir="results",
            databases=DatabasesConfig(
                ref_fa="ref.fa", genome="g.fa", gff3="g.gff3"
            ),
        )
        self._fm = FileManager(self._config)
        self._builder = RepIndexBuilder(self._config, self._fm)

    def test_returns_none_for_empty_transcripts(self):
        result = self._builder.select_representative("g1", [], {})
        assert result is None

    def test_returns_none_when_all_have_zero_cds(self):
        result = self._builder.select_representative("g1", ["t1", "t2"], {"t1": 0, "t2": 0})
        assert result is None

    def test_selects_longest_cds(self):
        result = self._builder.select_representative(
            "g1", ["t1", "t2"], {"t1": 100, "t2": 200}
        )
        assert result == "t2"

    def test_tie_broken_by_lexicographic_order(self):
        result = self._builder.select_representative(
            "g1", ["t1.2", "t1.1"], {"t1.1": 300, "t1.2": 300}
        )
        assert result == "t1.1"

    def test_ignores_transcripts_with_no_cds_entry(self):
        result = self._builder.select_representative(
            "g1", ["t1", "t2", "t3"], {"t1": 100, "t2": 50}
        )
        assert result == "t1"

    def test_single_transcript_selected(self):
        result = self._builder.select_representative("g1", ["t1"], {"t1": 150})
        assert result == "t1"


# ---------------------------------------------------------------------------
# run (integration)
# ---------------------------------------------------------------------------

class TestRepIndexBuilderRun:
    def test_run_creates_output_files(self, tmp_path):
        config = _make_config(tmp_path)
        fm = FileManager(config)
        builder = RepIndexBuilder(config, fm)
        builder.run()
        rep_dir = Path(config.genome_db.rep_index_dir)
        genome_name = config.project_name
        assert (rep_dir / f"{genome_name}.rep.cds.fa").exists()
        assert (rep_dir / f"{genome_name}.rep.pep.fa").exists()
        assert (rep_dir / f"{genome_name}.rep.idlist").exists()
        assert (rep_dir / f"{genome_name}.gene2rep.tsv").exists()

    def test_run_uses_genome_name(self, tmp_path):
        config = _make_config(tmp_path, genome_name="MyGenome")
        fm = FileManager(config)
        builder = RepIndexBuilder(config, fm)
        builder.run()
        rep_dir = Path(config.genome_db.rep_index_dir)
        assert (rep_dir / "MyGenome.rep.cds.fa").exists()

    def test_run_raises_when_gff3_missing(self, tmp_path):
        config = _make_config(tmp_path)
        Path(config.databases.gff3).unlink()
        fm = FileManager(config)
        builder = RepIndexBuilder(config, fm)
        with pytest.raises(StageInputError, match="GFF3"):
            builder.run()

    def test_run_raises_when_genome_missing(self, tmp_path):
        config = _make_config(tmp_path)
        Path(config.databases.genome).unlink()
        fm = FileManager(config)
        builder = RepIndexBuilder(config, fm)
        with pytest.raises(StageInputError, match="genome"):
            builder.run()

    def test_gene2rep_has_one_entry_per_gene(self, tmp_path):
        """Each gene must appear exactly once in gene2rep.tsv."""
        config = _make_config(tmp_path)
        fm = FileManager(config)
        builder = RepIndexBuilder(config, fm)
        builder.run()
        tsv = Path(config.genome_db.rep_index_dir) / f"{config.project_name}.gene2rep.tsv"
        lines = [l for l in tsv.read_text().splitlines() if l and not l.startswith("gene_id")]
        gene_ids = [l.split("\t")[0] for l in lines]
        assert len(gene_ids) == len(set(gene_ids))

    def test_rep_cds_contains_all_rep_transcripts(self, tmp_path):
        """All rep_transcript_ids in gene2rep.tsv must have sequences in rep.cds.fa."""
        config = _make_config(tmp_path)
        fm = FileManager(config)
        builder = RepIndexBuilder(config, fm)
        builder.run()
        rep_dir = Path(config.genome_db.rep_index_dir)
        genome_name = config.project_name

        tsv = rep_dir / f"{genome_name}.gene2rep.tsv"
        cds_fa = rep_dir / f"{genome_name}.rep.cds.fa"

        lines = [l for l in tsv.read_text().splitlines() if l and not l.startswith("gene_id")]
        rep_ids = {l.split("\t")[1] for l in lines}

        cds_content = cds_fa.read_text()
        for rep_id in rep_ids:
            assert f">{rep_id}" in cds_content, f"Missing CDS for {rep_id}"

    def test_selects_longest_cds_transcript(self, tmp_path):
        """gene1 has t1.1 (CDS 9bp) and t1.2 (CDS 15bp); t1.2 should be selected."""
        config = _make_config(tmp_path)
        fm = FileManager(config)
        builder = RepIndexBuilder(config, fm)
        builder.run()
        tsv = Path(config.genome_db.rep_index_dir) / f"{config.project_name}.gene2rep.tsv"
        lines = {l.split("\t")[0]: l.split("\t")[1]
                 for l in tsv.read_text().splitlines()
                 if l and not l.startswith("gene_id")}
        assert lines.get("gene1") == "t1.2"

    def test_gene_without_cds_skipped(self, tmp_path):
        """gene3 has no CDS; it should not appear in gene2rep.tsv."""
        config = _make_config(tmp_path)
        fm = FileManager(config)
        builder = RepIndexBuilder(config, fm)
        builder.run()
        tsv = Path(config.genome_db.rep_index_dir) / f"{config.project_name}.gene2rep.tsv"
        content = tsv.read_text()
        assert "gene3" not in content

    def test_run_skips_when_outputs_exist(self, tmp_path):
        config = _make_config(tmp_path)
        rep_dir = Path(config.genome_db.rep_index_dir)
        rep_dir.mkdir(parents=True, exist_ok=True)
        genome_name = config.project_name
        (rep_dir / f"{genome_name}.rep.cds.fa").write_text("existing")
        (rep_dir / f"{genome_name}.gene2rep.tsv").write_text("existing")
        fm = FileManager(config)
        builder = RepIndexBuilder(config, fm)
        builder.run(force=False)
        assert (rep_dir / f"{genome_name}.rep.cds.fa").read_text() == "existing"

    def test_run_force_overwrites(self, tmp_path):
        config = _make_config(tmp_path)
        rep_dir = Path(config.genome_db.rep_index_dir)
        rep_dir.mkdir(parents=True, exist_ok=True)
        genome_name = config.project_name
        (rep_dir / f"{genome_name}.rep.cds.fa").write_text("existing")
        (rep_dir / f"{genome_name}.gene2rep.tsv").write_text("existing")
        fm = FileManager(config)
        builder = RepIndexBuilder(config, fm)
        builder.run(force=True)
        assert (rep_dir / f"{genome_name}.rep.cds.fa").read_text() != "existing"

    def test_idlist_contains_rep_ids(self, tmp_path):
        config = _make_config(tmp_path)
        fm = FileManager(config)
        builder = RepIndexBuilder(config, fm)
        builder.run()
        rep_dir = Path(config.genome_db.rep_index_dir)
        genome_name = config.project_name
        idlist = (rep_dir / f"{genome_name}.rep.idlist").read_text().splitlines()
        tsv = rep_dir / f"{genome_name}.gene2rep.tsv"
        lines = [l for l in tsv.read_text().splitlines() if l and not l.startswith("gene_id")]
        rep_ids_from_tsv = {l.split("\t")[1] for l in lines}
        assert set(idlist) == rep_ids_from_tsv

    def test_gene2rep_has_correct_columns(self, tmp_path):
        config = _make_config(tmp_path)
        fm = FileManager(config)
        builder = RepIndexBuilder(config, fm)
        builder.run()
        tsv = Path(config.genome_db.rep_index_dir) / f"{config.project_name}.gene2rep.tsv"
        header = tsv.read_text().splitlines()[0]
        for col in ["gene_id", "rep_transcript_id", "cds_length", "transcript_count"]:
            assert col in header
