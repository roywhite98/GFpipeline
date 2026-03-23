"""Tests for CollinearityStage."""

from __future__ import annotations

import textwrap
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from gfpipeline.config.schema import (
    CollinearityConfig,
    DatabasesConfig,
    GenomeDbConfig,
    PipelineConfig,
)
from gfpipeline.core.exceptions import StageInputError
from gfpipeline.core.file_manager import FileManager
from gfpipeline.core.runner import ToolRunner
from gfpipeline.stages.collinearity import CollinearityStage


# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------

def _make_config(
    tmp_path: Path,
    tool: str = "jcvi",
    genome_name: str | None = None,
) -> PipelineConfig:
    return PipelineConfig(
        project_name="TEST",
        data_dir=str(tmp_path / "data"),
        result_dir=str(tmp_path / "results"),
        databases=DatabasesConfig(
            ref_fa=str(tmp_path / "ref.fa"),
            genome=str(tmp_path / "genome.fa"),
            gff3=str(tmp_path / "anno.gff3"),
        ),
        collinearity=CollinearityConfig(tool=tool, blast_threads=4),
        genome_db=GenomeDbConfig(genome_name=genome_name, rep_index_dir=str(tmp_path / "rep_index")),
    )


def _make_stage(
    tmp_path: Path,
    tool: str = "jcvi",
    genome_name: str | None = None,
    runner: ToolRunner | None = None,
) -> CollinearityStage:
    config = _make_config(tmp_path, tool=tool, genome_name=genome_name)
    fm = FileManager(config)
    fm.ensure_dirs()
    if runner is None:
        runner = MagicMock(spec=ToolRunner)
    return CollinearityStage(config, runner, fm)


_SAMPLE_GFF3 = textwrap.dedent("""\
    ##gff-version 3
    Chr1\t.\tgene\t1000\t2000\t.\t+\t.\tID=gene1;Name=gene1
    Chr1\t.\tmRNA\t1000\t2000\t.\t+\t.\tID=mRNA1;Parent=gene1
    Chr1\t.\tgene\t3000\t4000\t.\t-\t.\tID=gene2;Name=gene2
    Chr2\t.\tgene\t500\t1500\t.\t+\t.\tID=gene3;Name=gene3
    Chr2\t.\tgene\t2000\t3000\t.\t+\t.\tID=gene4;Name=gene4
""")

_SAMPLE_COLLINEARITY = textwrap.dedent("""\
    ## Alignment 0: score=150, e_value=1e-10, N=3, Chr1&Chr1 plus
    0-0\tgene1\tgene2\t1e-10
    0-1\tgene3\tgene4\t2e-10
    ## Alignment 1: score=80, e_value=1e-5, N=2, Chr2&Chr2 plus
    1-0\tgeneX\tgeneY\t1e-5
""")


# ---------------------------------------------------------------------------
# 1. parse_gene_locations with sample GFF3 content
# ---------------------------------------------------------------------------

def test_parse_gene_locations_basic(tmp_path):
    stage = _make_stage(tmp_path)
    gff3 = tmp_path / "anno.gff3"
    gff3.write_text(_SAMPLE_GFF3)

    stage.parse_gene_locations(gff3, ["gene1", "gene2", "gene3"])

    tsv = stage._gene_location_tsv
    assert tsv.exists()
    lines = tsv.read_text().splitlines()
    assert lines[0] == "gene_id\tchromosome\tstart\tend\tstrand"
    gene_ids_in_output = {line.split("\t")[0] for line in lines[1:] if line}
    assert gene_ids_in_output == {"gene1", "gene2", "gene3"}


# ---------------------------------------------------------------------------
# 2. extract_target_blocks with empty collinearity dir (writes header-only TSV)
# ---------------------------------------------------------------------------

def test_extract_target_blocks_empty_dir(tmp_path):
    stage = _make_stage(tmp_path)
    # collinearity dir does not exist → no files
    blocks = stage.extract_target_blocks(["gene1", "gene2"])

    assert blocks == []
    tsv = stage._blocks_tsv
    assert tsv.exists()
    lines = tsv.read_text().splitlines()
    assert lines[0] == "block_id\tgene1\tgene2\tchromosome1\tchromosome2\tscore"
    # Only header line
    assert len([l for l in lines if l.strip()]) == 1


# ---------------------------------------------------------------------------
# 3. run_all_vs_all_blast calls runner with correct command
# ---------------------------------------------------------------------------

def test_run_all_vs_all_blast_command(tmp_path):
    mock_runner = MagicMock(spec=ToolRunner)
    stage = _make_stage(tmp_path, runner=mock_runner)

    stage.run_all_vs_all_blast()

    mock_runner.run.assert_called_once()
    cmd = mock_runner.run.call_args[0][0]

    assert cmd[0] == "blastp"
    assert "-query" in cmd
    assert str(stage.fm.rep_pep) in cmd
    assert "-db" in cmd
    assert stage.fm.blast_db_prefix in cmd
    assert "-outfmt" in cmd
    assert "6" in cmd
    assert "-num_threads" in cmd
    assert "4" in cmd
    assert "-out" in cmd
    assert str(stage._blast_out) in cmd


# ---------------------------------------------------------------------------
# 4. run_jcvi calls runner with correct command
# ---------------------------------------------------------------------------

def test_run_jcvi_command(tmp_path):
    import sys
    mock_runner = MagicMock(spec=ToolRunner)
    stage = _make_stage(tmp_path, tool="jcvi", genome_name="rice", runner=mock_runner)

    stage.run_jcvi()

    mock_runner.run.assert_called_once()
    cmd = mock_runner.run.call_args[0][0]

    assert sys.executable in cmd
    assert "-m" in cmd
    assert "jcvi.compara.catalog" in cmd
    assert "ortholog" in cmd
    assert "rice" in cmd
    assert "--cpu=4" in cmd


# ---------------------------------------------------------------------------
# 5. run_mcscanx calls runner with correct command
# ---------------------------------------------------------------------------

def test_run_mcscanx_command(tmp_path):
    mock_runner = MagicMock(spec=ToolRunner)
    stage = _make_stage(tmp_path, tool="mcscanx", genome_name="rice", runner=mock_runner)

    stage.run_mcscanx()

    mock_runner.run.assert_called_once()
    cmd = mock_runner.run.call_args[0][0]

    assert cmd[0] == "MCScanX"
    # Input prefix should contain genome_name
    assert any("rice" in str(arg) for arg in cmd)


# ---------------------------------------------------------------------------
# 6. run raises StageInputError when gene idlist missing
# ---------------------------------------------------------------------------

def test_run_raises_if_gene_idlist_missing(tmp_path):
    stage = _make_stage(tmp_path)
    # Create gff3 and rep_pep but NOT the gene idlist
    Path(stage.config.databases.gff3).write_text(_SAMPLE_GFF3)
    stage.fm.rep_pep.parent.mkdir(parents=True, exist_ok=True)
    stage.fm.rep_pep.write_text(">gene1\nMACDEFGH\n")

    with pytest.raises(StageInputError):
        stage.run()


# ---------------------------------------------------------------------------
# 7. run raises StageInputError when gff3 missing
# ---------------------------------------------------------------------------

def test_run_raises_if_gff3_missing(tmp_path):
    stage = _make_stage(tmp_path)
    # Create gene idlist and rep_pep but NOT gff3
    stage.fm.ensure_dirs()
    stage._gene_idlist.write_text("gene1\n")
    stage.fm.rep_pep.parent.mkdir(parents=True, exist_ok=True)
    stage.fm.rep_pep.write_text(">gene1\nMACDEFGH\n")

    with pytest.raises(StageInputError):
        stage.run()


# ---------------------------------------------------------------------------
# 8. run raises StageInputError when pep missing
# ---------------------------------------------------------------------------

def test_run_raises_if_pep_missing(tmp_path):
    stage = _make_stage(tmp_path)
    # Create gene idlist and gff3 but NOT rep_pep
    stage.fm.ensure_dirs()
    stage._gene_idlist.write_text("gene1\n")
    Path(stage.config.databases.gff3).write_text(_SAMPLE_GFF3)

    with pytest.raises(StageInputError):
        stage.run()


# ---------------------------------------------------------------------------
# 9. run skips if outputs exist and force=False
# ---------------------------------------------------------------------------

def test_run_skips_if_outputs_exist(tmp_path):
    mock_runner = MagicMock(spec=ToolRunner)
    stage = _make_stage(tmp_path, runner=mock_runner)

    # Create all required inputs
    stage.fm.ensure_dirs()
    stage._gene_idlist.write_text("gene1\n")
    Path(stage.config.databases.gff3).write_text(_SAMPLE_GFF3)
    stage.fm.rep_pep.parent.mkdir(parents=True, exist_ok=True)
    stage.fm.rep_pep.write_text(">gene1\nMACDEFGH\n")

    # Create both outputs
    stage._blocks_tsv.write_text("block_id\tgene1\tgene2\tchromosome1\tchromosome2\tscore\n")
    stage._gene_location_tsv.write_text("gene_id\tchromosome\tstart\tend\tstrand\n")

    stage.run(force=False)

    mock_runner.run.assert_not_called()


# ---------------------------------------------------------------------------
# 10. parse_gene_locations only includes genes in the provided list
# ---------------------------------------------------------------------------

def test_parse_gene_locations_filters_by_list(tmp_path):
    stage = _make_stage(tmp_path)
    gff3 = tmp_path / "anno.gff3"
    gff3.write_text(_SAMPLE_GFF3)

    # Only request gene1 and gene4
    stage.parse_gene_locations(gff3, ["gene1", "gene4"])

    tsv = stage._gene_location_tsv
    lines = tsv.read_text().splitlines()
    gene_ids_in_output = {line.split("\t")[0] for line in lines[1:] if line}
    assert gene_ids_in_output == {"gene1", "gene4"}
    # gene2 and gene3 must NOT be present
    assert "gene2" not in gene_ids_in_output
    assert "gene3" not in gene_ids_in_output
