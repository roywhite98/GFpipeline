"""Unit tests for stages/tree.py — TreeStage."""

from __future__ import annotations

from pathlib import Path

import pytest

from gfpipeline.config.schema import DatabasesConfig, PipelineConfig, TreeConfig
from gfpipeline.core.exceptions import StageInputError
from gfpipeline.core.file_manager import FileManager
from gfpipeline.core.runner import ToolRunner
from gfpipeline.stages.tree import TreeStage


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def make_config(tmp_path: Path, bootstrap: int = 1000) -> PipelineConfig:
    data_dir = tmp_path / "data"
    result_dir = tmp_path / "results"
    data_dir.mkdir()
    result_dir.mkdir()
    return PipelineConfig(
        project_name="TEST",
        data_dir=str(data_dir),
        result_dir=str(result_dir),
        databases=DatabasesConfig(
            ref_fa=str(tmp_path / "ref.fa"),
            genome=str(tmp_path / "genome.fa"),
            gff3=str(tmp_path / "anno.gff3"),
        ),
        tree=TreeConfig(iqtree_bootstrap=bootstrap),
    )


def make_stage(tmp_path: Path, bootstrap: int = 1000) -> TreeStage:
    config = make_config(tmp_path, bootstrap)
    runner = ToolRunner(dry_run=True)
    fm = FileManager(config)
    return TreeStage(config, runner, fm)


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

class TestRunInputValidation:
    def test_raises_when_pep_fa_missing(self, tmp_path):
        stage = make_stage(tmp_path)
        with pytest.raises(StageInputError) as exc_info:
            stage.run()
        assert "identify" in str(exc_info.value).lower() or "candidates.pep" in str(exc_info.value)

    def test_no_error_when_pep_fa_exists(self, tmp_path):
        stage = make_stage(tmp_path)
        # Create the required input file
        stage._pep_fa.write_text(">seq1\nMAAA\n")
        # dry_run=True so no real tools are called
        stage.run()


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

class TestPathHelpers:
    def test_pep_fa_path(self, tmp_path):
        stage = make_stage(tmp_path)
        expected = Path(stage.fm.result_dir) / "TEST.identify.candidates.pep.fa"
        assert stage._pep_fa == expected

    def test_tree_afa_path(self, tmp_path):
        stage = make_stage(tmp_path)
        expected = Path(stage.fm.result_dir) / "TEST.tree.pep.afa"
        assert stage._tree_afa == expected

    def test_tree_trimed_afa_path(self, tmp_path):
        stage = make_stage(tmp_path)
        expected = Path(stage.fm.result_dir) / "TEST.tree.pep.trimed.afa"
        assert stage._tree_trimed_afa == expected


# ---------------------------------------------------------------------------
# Tool invocation order (dry_run captures commands)
# ---------------------------------------------------------------------------

class TestToolInvocationOrder:
    def test_muscle_trimal_iqtree_called_in_order(self, tmp_path):
        stage = make_stage(tmp_path)
        stage._pep_fa.write_text(">seq1\nMAAA\n")

        called_tools = []
        original_run = stage.runner.run

        def capture_run(cmd, **kwargs):
            called_tools.append(cmd[0])
            return original_run(cmd, **kwargs)

        stage.runner.run = capture_run
        stage.run()

        assert called_tools[0] == "muscle"
        assert called_tools[1] == "trimal"
        assert called_tools[2] == "iqtree2"

    def test_muscle_skipped_if_afa_exists(self, tmp_path):
        stage = make_stage(tmp_path)
        stage._pep_fa.write_text(">seq1\nMAAA\n")
        # Pre-create the afa output
        stage._tree_afa.write_text("dummy afa")

        called_tools = []
        original_run = stage.runner.run

        def capture_run(cmd, **kwargs):
            called_tools.append(cmd[0])
            return original_run(cmd, **kwargs)

        stage.runner.run = capture_run
        stage.run()

        assert "muscle" not in called_tools
        assert "trimal" in called_tools
        assert "iqtree2" in called_tools

    def test_trimal_skipped_if_trimed_afa_exists(self, tmp_path):
        stage = make_stage(tmp_path)
        stage._pep_fa.write_text(">seq1\nMAAA\n")
        stage._tree_afa.write_text("dummy afa")
        stage._tree_trimed_afa.write_text("dummy trimed afa")

        called_tools = []
        original_run = stage.runner.run

        def capture_run(cmd, **kwargs):
            called_tools.append(cmd[0])
            return original_run(cmd, **kwargs)

        stage.runner.run = capture_run
        stage.run()

        assert "muscle" not in called_tools
        assert "trimal" not in called_tools
        assert "iqtree2" in called_tools


# ---------------------------------------------------------------------------
# IQ-TREE bootstrap parameter
# ---------------------------------------------------------------------------

class TestIqtreeBootstrap:
    def test_bootstrap_value_passed_to_iqtree(self, tmp_path):
        stage = make_stage(tmp_path, bootstrap=500)
        stage._pep_fa.write_text(">seq1\nMAAA\n")

        captured_cmds = []
        original_run = stage.runner.run

        def capture_run(cmd, **kwargs):
            captured_cmds.append(cmd)
            return original_run(cmd, **kwargs)

        stage.runner.run = capture_run
        stage.run()

        iqtree_cmd = next(c for c in captured_cmds if c[0] == "iqtree2")
        bb_idx = iqtree_cmd.index("-bb")
        assert iqtree_cmd[bb_idx + 1] == "500"

    def test_iqtree_uses_bnni_and_nt_auto(self, tmp_path):
        stage = make_stage(tmp_path)
        stage._pep_fa.write_text(">seq1\nMAAA\n")

        captured_cmds = []
        original_run = stage.runner.run

        def capture_run(cmd, **kwargs):
            captured_cmds.append(cmd)
            return original_run(cmd, **kwargs)

        stage.runner.run = capture_run
        stage.run()

        iqtree_cmd = next(c for c in captured_cmds if c[0] == "iqtree2")
        assert "-bnni" in iqtree_cmd
        assert "-nt" in iqtree_cmd
        assert "AUTO" in iqtree_cmd

    def test_iqtree_input_is_trimed_afa(self, tmp_path):
        stage = make_stage(tmp_path)
        stage._pep_fa.write_text(">seq1\nMAAA\n")

        captured_cmds = []
        original_run = stage.runner.run

        def capture_run(cmd, **kwargs):
            captured_cmds.append(cmd)
            return original_run(cmd, **kwargs)

        stage.runner.run = capture_run
        stage.run()

        iqtree_cmd = next(c for c in captured_cmds if c[0] == "iqtree2")
        s_idx = iqtree_cmd.index("-s")
        assert "trimed" in iqtree_cmd[s_idx + 1]
