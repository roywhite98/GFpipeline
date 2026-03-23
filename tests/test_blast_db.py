"""Unit tests for genome_db/blast_db.py — BlastDbBuilder."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from gfpipeline.config.schema import (
    DatabasesConfig,
    GenomeDbConfig,
    PipelineConfig,
    ToolsConfig,
)
from gfpipeline.core.exceptions import PipelineError, StageInputError
from gfpipeline.core.file_manager import FileManager
from gfpipeline.core.runner import ToolRunner
from gfpipeline.genome_db.blast_db import BlastDbBuilder, _db_exists


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_config(tmp_path: Path, build_prot: bool = True, build_nucl: bool = True) -> PipelineConfig:
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nATCG\n")

    # Create rep_pep file in the expected rep_index_dir location
    rep_dir = tmp_path / "rep_index"
    rep_dir.mkdir(parents=True, exist_ok=True)
    (rep_dir / "TEST.rep.pep.fa").write_text(">seq1\nMACK\n")

    return PipelineConfig(
        project_name="TEST",
        data_dir=str(tmp_path / "data"),
        result_dir=str(tmp_path / "results"),
        databases=DatabasesConfig(
            ref_fa=str(tmp_path / "ref.fa"),
            genome=str(genome),
            gff3=str(tmp_path / "genome.gff3"),
        ),
        genome_db=GenomeDbConfig(
            build_prot=build_prot,
            build_nucl=build_nucl,
            rep_index_dir=str(rep_dir),
        ),
    )


@pytest.fixture
def dry_runner():
    return ToolRunner(dry_run=True)


# ---------------------------------------------------------------------------
# _db_exists helper
# ---------------------------------------------------------------------------

class TestDbExists:
    def test_returns_false_when_files_missing(self, tmp_path):
        assert _db_exists(str(tmp_path / "db"), [".phr", ".pin"]) is False

    def test_returns_true_when_all_files_present(self, tmp_path):
        prefix = str(tmp_path / "db")
        for ext in [".phr", ".pin"]:
            Path(prefix + ext).write_text("")
        assert _db_exists(prefix, [".phr", ".pin"]) is True

    def test_returns_false_when_only_some_files_present(self, tmp_path):
        prefix = str(tmp_path / "db")
        Path(prefix + ".phr").write_text("")
        assert _db_exists(prefix, [".phr", ".pin"]) is False


# ---------------------------------------------------------------------------
# BlastDbBuilder.run
# ---------------------------------------------------------------------------

class TestBlastDbBuilderRun:
    def test_raises_when_both_disabled(self, tmp_path):
        config = _make_config(tmp_path, build_prot=False, build_nucl=False)
        runner = ToolRunner(dry_run=True)
        fm = FileManager(config)
        builder = BlastDbBuilder(config, runner, fm)
        with pytest.raises(PipelineError, match="build_prot and build_nucl are False"):
            builder.run()

    def test_raises_when_pep_missing(self, tmp_path):
        config = _make_config(tmp_path, build_prot=True, build_nucl=False)
        # Remove rep_pep file
        fm = FileManager(config)
        fm.rep_pep.unlink()
        runner = ToolRunner(dry_run=True)
        builder = BlastDbBuilder(config, runner, fm)
        with pytest.raises(StageInputError, match="protein FASTA not found"):
            builder.run()

    def test_raises_when_genome_missing(self, tmp_path):
        config = _make_config(tmp_path, build_prot=False, build_nucl=True)
        Path(config.databases.genome).unlink()
        runner = ToolRunner(dry_run=True)
        fm = FileManager(config)
        builder = BlastDbBuilder(config, runner, fm)
        with pytest.raises(StageInputError, match="genome FASTA not found"):
            builder.run()

    def test_dry_run_does_not_fail(self, tmp_path):
        config = _make_config(tmp_path)
        runner = ToolRunner(dry_run=True)
        fm = FileManager(config)
        builder = BlastDbBuilder(config, runner, fm)
        # Should not raise; dry_run skips actual makeblastdb
        builder.run()

    def test_skips_prot_db_when_exists(self, tmp_path):
        config = _make_config(tmp_path, build_prot=True, build_nucl=False)
        fm = FileManager(config)
        prefix = fm.blast_db_prefix
        for ext in [".phr", ".pin", ".psq"]:
            Path(prefix + ext).write_text("")

        mock_runner = MagicMock()
        builder = BlastDbBuilder(config, mock_runner, fm)
        builder.run(force=False)
        mock_runner.run.assert_not_called()

    def test_force_rebuilds_existing_db(self, tmp_path):
        config = _make_config(tmp_path, build_prot=True, build_nucl=False)
        fm = FileManager(config)
        prefix = fm.blast_db_prefix
        for ext in [".phr", ".pin", ".psq"]:
            Path(prefix + ext).write_text("")

        runner = ToolRunner(dry_run=True)
        builder = BlastDbBuilder(config, runner, fm)
        # Should not raise even with force=True (dry_run)
        builder.run(force=True)

    def test_calls_makeblastdb_for_prot(self, tmp_path):
        config = _make_config(tmp_path, build_prot=True, build_nucl=False)
        mock_runner = MagicMock()
        fm = FileManager(config)
        builder = BlastDbBuilder(config, mock_runner, fm)
        builder.run()
        mock_runner.run.assert_called_once()
        call_args = mock_runner.run.call_args[0][0]
        assert "makeblastdb" in call_args[0]
        assert "-dbtype" in call_args
        assert "prot" in call_args

    def test_calls_makeblastdb_for_nucl(self, tmp_path):
        config = _make_config(tmp_path, build_prot=False, build_nucl=True)
        mock_runner = MagicMock()
        fm = FileManager(config)
        builder = BlastDbBuilder(config, mock_runner, fm)
        builder.run()
        mock_runner.run.assert_called_once()
        call_args = mock_runner.run.call_args[0][0]
        assert "nucl" in call_args

    def test_calls_makeblastdb_for_both(self, tmp_path):
        config = _make_config(tmp_path, build_prot=True, build_nucl=True)
        mock_runner = MagicMock()
        fm = FileManager(config)
        builder = BlastDbBuilder(config, mock_runner, fm)
        builder.run()
        assert mock_runner.run.call_count == 2
