"""Tests for TransStage."""

from __future__ import annotations

import textwrap
from pathlib import Path

import pandas as pd
import pytest

from gfpipeline.config.schema import DatabasesConfig, PipelineConfig, TransConfig
from gfpipeline.core.exceptions import StageInputError
from gfpipeline.core.file_manager import FileManager
from gfpipeline.stages.trans import TransStage


# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------

def _make_config(tmp_path: Path, trans: TransConfig | None = None) -> PipelineConfig:
    return PipelineConfig(
        project_name="TEST",
        data_dir=str(tmp_path / "data"),
        result_dir=str(tmp_path / "results"),
        databases=DatabasesConfig(
            ref_fa=str(tmp_path / "ref.fa"),
            genome=str(tmp_path / "genome.fa"),
            gff3=str(tmp_path / "anno.gff3"),
        ),
        trans=trans or TransConfig(),
    )


def _make_stage(tmp_path: Path, trans: TransConfig | None = None) -> TransStage:
    config = _make_config(tmp_path, trans)
    fm = FileManager(config)
    fm.ensure_dirs()
    return TransStage(config, fm)


_TSV_CONTENT = textwrap.dedent("""\
    gene_id\tsample1\tsample2\tsample3
    GeneA\t1.0\t2.0\t3.0
    GeneB\t0.0\t4.0\t5.0
    GeneC\t6.0\t0.0\t7.0
""")


def _write_matrix(tmp_path: Path, content: str = _TSV_CONTENT) -> Path:
    p = tmp_path / "expr.tsv"
    p.write_text(content)
    return p


# ---------------------------------------------------------------------------
# 1. load_expression_matrix with a simple TSV file
# ---------------------------------------------------------------------------

def test_load_expression_matrix_tsv(tmp_path):
    stage = _make_stage(tmp_path)
    matrix_file = _write_matrix(tmp_path)

    df = stage.load_expression_matrix(matrix_file)

    assert isinstance(df, pd.DataFrame)
    assert list(df.index) == ["GeneA", "GeneB", "GeneC"]
    assert list(df.columns) == ["sample1", "sample2", "sample3"]
    assert df.loc["GeneA", "sample1"] == 1.0


# ---------------------------------------------------------------------------
# 2. filter_family_members returns only matching rows
# ---------------------------------------------------------------------------

def test_filter_family_members_returns_matching(tmp_path):
    stage = _make_stage(tmp_path)
    matrix_file = _write_matrix(tmp_path)
    df = stage.load_expression_matrix(matrix_file)

    result = stage.filter_family_members(df, ["GeneA", "GeneC"])

    assert list(result.index) == ["GeneA", "GeneC"]
    assert "GeneB" not in result.index


# ---------------------------------------------------------------------------
# 3. filter_family_members returns empty DataFrame when no matches
# ---------------------------------------------------------------------------

def test_filter_family_members_no_matches(tmp_path):
    stage = _make_stage(tmp_path)
    matrix_file = _write_matrix(tmp_path)
    df = stage.load_expression_matrix(matrix_file)

    result = stage.filter_family_members(df, ["GeneX", "GeneY"])

    assert result.empty


# ---------------------------------------------------------------------------
# 4. run raises StageInputError when expression_matrix is None
# ---------------------------------------------------------------------------

def test_run_raises_when_expression_matrix_none(tmp_path):
    stage = _make_stage(tmp_path, TransConfig(expression_matrix=None))

    with pytest.raises(StageInputError):
        stage.run()


# ---------------------------------------------------------------------------
# 5. run raises StageInputError when expression_matrix file doesn't exist
# ---------------------------------------------------------------------------

def test_run_raises_when_expression_matrix_missing(tmp_path):
    stage = _make_stage(
        tmp_path,
        TransConfig(expression_matrix=str(tmp_path / "nonexistent.tsv")),
    )

    with pytest.raises(StageInputError):
        stage.run()


# ---------------------------------------------------------------------------
# 6. run skips if outputs exist and force=False
# ---------------------------------------------------------------------------

def test_run_skips_if_outputs_exist(tmp_path):
    matrix_file = _write_matrix(tmp_path)
    stage = _make_stage(tmp_path, TransConfig(expression_matrix=str(matrix_file)))

    # Pre-create all output files
    for p in [stage._heatmap_pdf, stage._deg_tsv, stage._venn_pdf]:
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("existing\n")

    stage.run(force=False)

    # Outputs should be unchanged
    assert stage._deg_tsv.read_text() == "existing\n"


# ---------------------------------------------------------------------------
# 7. filter_deg with logfc_threshold set
# ---------------------------------------------------------------------------

def test_filter_deg_with_logfc_threshold(tmp_path):
    stage = _make_stage(tmp_path, TransConfig(logfc_threshold=3.0))
    matrix_file = _write_matrix(tmp_path)
    df = stage.load_expression_matrix(matrix_file)

    result = stage.filter_deg(df)

    # GeneA max=3.0, GeneB max=5.0, GeneC max=7.0 → all >= 3.0
    # GeneA has a value of exactly 3.0 which passes abs >= 3.0
    assert "GeneB" in result.index
    assert "GeneC" in result.index


def test_filter_deg_logfc_threshold_excludes_low(tmp_path):
    stage = _make_stage(tmp_path, TransConfig(logfc_threshold=5.0))
    matrix_file = _write_matrix(tmp_path)
    df = stage.load_expression_matrix(matrix_file)

    result = stage.filter_deg(df)

    # GeneA max=3.0 < 5.0 → excluded
    assert "GeneA" not in result.index
    # GeneB max=5.0 >= 5.0 → included
    assert "GeneB" in result.index
    # GeneC max=7.0 >= 5.0 → included
    assert "GeneC" in result.index


# ---------------------------------------------------------------------------
# 8. filter_deg returns full matrix when no thresholds set
# ---------------------------------------------------------------------------

def test_filter_deg_no_thresholds_returns_full(tmp_path):
    stage = _make_stage(tmp_path, TransConfig())
    matrix_file = _write_matrix(tmp_path)
    df = stage.load_expression_matrix(matrix_file)

    result = stage.filter_deg(df)

    assert len(result) == len(df)
    assert list(result.index) == list(df.index)


# ---------------------------------------------------------------------------
# 9. plot_heatmap creates a PDF file
# ---------------------------------------------------------------------------

def test_plot_heatmap_creates_pdf(tmp_path):
    stage = _make_stage(tmp_path)
    matrix_file = _write_matrix(tmp_path)
    df = stage.load_expression_matrix(matrix_file)

    output = tmp_path / "heatmap.pdf"
    stage.plot_heatmap(df, output)

    assert output.exists()
    assert output.stat().st_size > 0


# ---------------------------------------------------------------------------
# 10. run logs warning when no family members in matrix
# ---------------------------------------------------------------------------

def test_run_warns_when_no_family_members(tmp_path, caplog):
    import logging

    matrix_file = _write_matrix(tmp_path)
    # gene idlist with IDs not in the matrix
    stage = _make_stage(tmp_path, TransConfig(expression_matrix=str(matrix_file)))
    stage._gene_idlist.parent.mkdir(parents=True, exist_ok=True)
    stage._gene_idlist.write_text("GeneX\nGeneY\n")

    with caplog.at_level(logging.WARNING):
        stage.run(force=True)

    assert any("No gene family members" in r.message for r in caplog.records)
