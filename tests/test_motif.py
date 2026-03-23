"""Tests for MotifStage."""

from __future__ import annotations

import textwrap
from pathlib import Path
from unittest.mock import MagicMock, call

import pytest

from gfpipeline.config.schema import (
    DatabasesConfig,
    GenomeDbConfig,
    MotifConfig,
    PipelineConfig,
)
from gfpipeline.core.exceptions import StageInputError
from gfpipeline.core.file_manager import FileManager
from gfpipeline.core.runner import ToolRunner
from gfpipeline.stages.motif import MotifStage


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------

def _make_config(tmp_path: Path, motif_cfg: MotifConfig | None = None) -> PipelineConfig:
    return PipelineConfig(
        project_name="TEST",
        data_dir=str(tmp_path / "data"),
        result_dir=str(tmp_path / "results"),
        databases=DatabasesConfig(
            ref_fa=str(tmp_path / "ref.fa"),
            genome=str(tmp_path / "genome.fa"),
            gff3=str(tmp_path / "anno.gff3"),
        ),
        motif=motif_cfg or MotifConfig(),
        genome_db=GenomeDbConfig(rep_index_dir=str(tmp_path / "rep_index")),
    )


def _make_stage(
    tmp_path: Path,
    motif_cfg: MotifConfig | None = None,
    runner: ToolRunner | None = None,
) -> MotifStage:
    config = _make_config(tmp_path, motif_cfg)
    fm = FileManager(config)
    fm.ensure_dirs()
    if runner is None:
        runner = MagicMock(spec=ToolRunner)
    return MotifStage(config, runner, fm)


# Sample fimo.tsv content (tab-separated)
_FIMO_TSV = textwrap.dedent("""\
    #FIMO version 5.5.0
    motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence
    MEME-1\tmotif_A\tgene1\t10\t25\t+\t15.3\t1.2e-5\t0.01\tACDEFGHIKL
    MEME-2\tmotif_B\tgene1\t50\t65\t+\t12.1\t3.4e-4\t0.05\tMNPQRSTVWY
    MEME-1\tmotif_A\tgene2\t5\t20\t+\t14.0\t2.0e-5\t0.02\tACDEFGHIKL
    MEME-3\tmotif_C\tgene3\t1\t16\t+\t10.5\t5.0e-3\t0.1\tGHIKLMNPQR
""")


# ---------------------------------------------------------------------------
# 1. parse_fimo with sample content
# ---------------------------------------------------------------------------

def test_parse_fimo_basic(tmp_path):
    stage = _make_stage(tmp_path)
    fimo_dir = tmp_path / "fimo_out"
    fimo_dir.mkdir()
    (fimo_dir / "fimo.tsv").write_text(_FIMO_TSV)

    result = stage.parse_fimo(fimo_dir)

    assert "gene1" in result
    assert "gene2" in result
    assert "gene3" in result
    assert "MEME-1" in result["gene1"]
    assert "MEME-2" in result["gene1"]
    assert result["gene2"] == ["MEME-1"]
    assert result["gene3"] == ["MEME-3"]


# ---------------------------------------------------------------------------
# 2. filter_genes mode=any
# ---------------------------------------------------------------------------

def test_filter_genes_any(tmp_path):
    stage = _make_stage(tmp_path, MotifConfig(filter_mode="any"))
    gene_motifs = {
        "gene1": ["MEME-1", "MEME-2"],
        "gene2": ["MEME-1"],
        "gene3": ["MEME-3"],
    }
    result = stage.filter_genes(gene_motifs, ["MEME-1", "MEME-2", "MEME-3"])
    assert sorted(result) == ["gene1", "gene2", "gene3"]


# ---------------------------------------------------------------------------
# 3. filter_genes mode=all
# ---------------------------------------------------------------------------

def test_filter_genes_all(tmp_path):
    stage = _make_stage(tmp_path, MotifConfig(filter_mode="all"))
    gene_motifs = {
        "gene1": ["MEME-1", "MEME-2"],
        "gene2": ["MEME-1"],
        "gene3": ["MEME-3"],
    }
    result = stage.filter_genes(gene_motifs, ["MEME-1", "MEME-2"])
    # Only gene1 has both MEME-1 and MEME-2
    assert result == ["gene1"]


# ---------------------------------------------------------------------------
# 4. filter_genes mode=min_count
# ---------------------------------------------------------------------------

def test_filter_genes_min_count(tmp_path):
    stage = _make_stage(tmp_path, MotifConfig(filter_mode="min_count", min_motif_count=2))
    gene_motifs = {
        "gene1": ["MEME-1", "MEME-2"],
        "gene2": ["MEME-1"],
        "gene3": ["MEME-1", "MEME-2", "MEME-3"],
    }
    result = stage.filter_genes(gene_motifs, ["MEME-1", "MEME-2", "MEME-3"])
    # gene1 has 2, gene3 has 3 — both >= 2
    assert sorted(result) == ["gene1", "gene3"]


# ---------------------------------------------------------------------------
# 5. filter_genes with empty gene_motifs
# ---------------------------------------------------------------------------

def test_filter_genes_empty(tmp_path):
    stage = _make_stage(tmp_path, MotifConfig(filter_mode="any"))
    result = stage.filter_genes({}, [])
    assert result == []


def test_filter_genes_all_empty(tmp_path):
    stage = _make_stage(tmp_path, MotifConfig(filter_mode="all"))
    result = stage.filter_genes({}, ["MEME-1"])
    assert result == []


def test_filter_genes_min_count_empty(tmp_path):
    stage = _make_stage(tmp_path, MotifConfig(filter_mode="min_count", min_motif_count=1))
    result = stage.filter_genes({}, [])
    assert result == []


# ---------------------------------------------------------------------------
# 6. run raises StageInputError when candidates.pep.fa missing
# ---------------------------------------------------------------------------

def test_run_raises_if_pep_missing(tmp_path):
    stage = _make_stage(tmp_path)
    with pytest.raises(StageInputError):
        stage.run()


# ---------------------------------------------------------------------------
# 7. run skips if outputs exist and force=False
# ---------------------------------------------------------------------------

def test_run_skips_if_outputs_exist(tmp_path):
    mock_runner = MagicMock(spec=ToolRunner)
    stage = _make_stage(tmp_path, runner=mock_runner)

    # Create required input and both outputs
    stage._pep_fa.write_text(">gene1\nMACDEFGH\n")
    stage._candidates_idlist.write_text("gene1\n")
    stage._summary_tsv.write_text("gene_id\tmotif_id\n")

    stage.run(force=False)

    mock_runner.run.assert_not_called()


# ---------------------------------------------------------------------------
# 8. run_meme calls runner with correct command
# ---------------------------------------------------------------------------

def test_run_meme_correct_command(tmp_path):
    mock_runner = MagicMock(spec=ToolRunner)
    motif_cfg = MotifConfig(num_motifs=5, min_width=8, max_width=40)
    stage = _make_stage(tmp_path, motif_cfg, runner=mock_runner)

    # Create the pep file so the path exists
    stage._pep_fa.write_text(">gene1\nMACDEFGH\n")

    # run_meme will call runner.run; we just check the command
    stage.run_meme()

    mock_runner.run.assert_called_once()
    cmd = mock_runner.run.call_args[0][0]

    assert cmd[0] == "meme"
    assert str(stage._pep_fa) in cmd
    assert "-protein" in cmd
    assert "-nmotifs" in cmd
    assert "5" in cmd
    assert "-minw" in cmd
    assert "8" in cmd
    assert "-maxw" in cmd
    assert "40" in cmd
    assert "-oc" in cmd
    assert str(stage._meme_dir) in cmd


# ---------------------------------------------------------------------------
# 9. run_fimo calls runner with correct command
# ---------------------------------------------------------------------------

def test_run_fimo_correct_command(tmp_path):
    mock_runner = MagicMock(spec=ToolRunner)
    stage = _make_stage(tmp_path, runner=mock_runner)

    meme_dir = tmp_path / "results" / "TEST.motif.meme"
    meme_dir.mkdir(parents=True)
    (meme_dir / "meme.xml").write_text("<MEME/>")

    # Create rep_pep db
    stage.fm.rep_pep.parent.mkdir(parents=True, exist_ok=True)
    stage.fm.rep_pep.write_text(">g1\nMACDEFGH\n")

    stage.run_fimo(meme_dir)

    mock_runner.run.assert_called_once()
    cmd = mock_runner.run.call_args[0][0]

    assert cmd[0] == "fimo"
    assert "--oc" in cmd
    assert str(stage._fimo_dir) in cmd
    assert str(meme_dir / "meme.xml") in cmd
    assert str(stage.fm.rep_pep) in cmd


# ---------------------------------------------------------------------------
# 10. parse_fimo handles comment lines and header lines correctly
# ---------------------------------------------------------------------------

def test_parse_fimo_skips_comments_and_header(tmp_path):
    stage = _make_stage(tmp_path)
    fimo_dir = tmp_path / "fimo_out"
    fimo_dir.mkdir()

    content = textwrap.dedent("""\
        # This is a comment
        # Another comment
        motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence
        MEME-1\tmotif_A\tgeneX\t1\t10\t+\t9.0\t1e-4\t0.05\tACDEFGHIKL
        # trailing comment
    """)
    (fimo_dir / "fimo.tsv").write_text(content)

    result = stage.parse_fimo(fimo_dir)

    # Only geneX should be present; header and comments must be skipped
    assert list(result.keys()) == ["geneX"]
    assert result["geneX"] == ["MEME-1"]


def test_parse_fimo_no_duplicate_motifs_per_gene(tmp_path):
    """Same motif appearing multiple times for a gene should be deduplicated."""
    stage = _make_stage(tmp_path)
    fimo_dir = tmp_path / "fimo_out"
    fimo_dir.mkdir()

    content = textwrap.dedent("""\
        motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence
        MEME-1\tmotif_A\tgene1\t1\t10\t+\t9.0\t1e-4\t0.05\tACDEFGHIKL
        MEME-1\tmotif_A\tgene1\t20\t30\t+\t8.5\t2e-4\t0.06\tACDEFGHIKL
    """)
    (fimo_dir / "fimo.tsv").write_text(content)

    result = stage.parse_fimo(fimo_dir)
    assert result["gene1"].count("MEME-1") == 1


def test_filter_genes_returns_sorted(tmp_path):
    """filter_genes should return a sorted list."""
    stage = _make_stage(tmp_path, MotifConfig(filter_mode="any"))
    gene_motifs = {"geneZ": ["MEME-1"], "geneA": ["MEME-1"], "geneM": ["MEME-1"]}
    result = stage.filter_genes(gene_motifs, ["MEME-1"])
    assert result == sorted(result)
