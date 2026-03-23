"""Unit tests for core/file_manager.py."""

from pathlib import Path

import pytest

from gfpipeline.config.schema import DatabasesConfig, PipelineConfig
from gfpipeline.core.file_manager import FileManager


def make_config(tmp_path: Path, proj: str = "OSACO") -> PipelineConfig:
    return PipelineConfig(
        project_name=proj,
        data_dir=str(tmp_path / "data"),
        result_dir=str(tmp_path / "results"),
        databases=DatabasesConfig(
            ref_fa="ref.fa",
            genome="genome.fa",
            gff3="anno.gff3",
        ),
    )


class TestFileManagerResult:
    def test_result_path_format(self, tmp_path):
        fm = FileManager(make_config(tmp_path, "PROJ"))
        p = fm.result("identify", "hmm", "out")
        assert p.name == "PROJ.identify.hmm.out"

    def test_result_path_is_under_result_dir(self, tmp_path):
        fm = FileManager(make_config(tmp_path))
        p = fm.result("tree", "pep", "afa")
        assert p.parent == Path(str(tmp_path / "results"))

    def test_result_path_components(self, tmp_path):
        fm = FileManager(make_config(tmp_path, "MyProj"))
        p = fm.result("domain", "cdd", "txt")
        assert p.name == "MyProj.domain.cdd.txt"


class TestFileManagerData:
    def test_data_path_is_under_data_dir(self, tmp_path):
        fm = FileManager(make_config(tmp_path))
        p = fm.data("OSACO.hmm")
        assert p == Path(str(tmp_path / "data")) / "OSACO.hmm"

    def test_data_path_preserves_name(self, tmp_path):
        fm = FileManager(make_config(tmp_path))
        assert fm.data("ref.fa").name == "ref.fa"


class TestFileManagerEnsureDirs:
    def test_ensure_dirs_creates_data_dir(self, tmp_path):
        fm = FileManager(make_config(tmp_path))
        fm.ensure_dirs()
        assert (tmp_path / "data").exists()

    def test_ensure_dirs_creates_result_dir(self, tmp_path):
        fm = FileManager(make_config(tmp_path))
        fm.ensure_dirs()
        assert (tmp_path / "results").exists()

    def test_ensure_dirs_idempotent(self, tmp_path):
        fm = FileManager(make_config(tmp_path))
        fm.ensure_dirs()
        fm.ensure_dirs()  # should not raise


class TestFileManagerSkipIfExists:
    def test_skip_if_exists_returns_true_when_file_exists_no_force(self, tmp_path):
        fm = FileManager(make_config(tmp_path))
        f = tmp_path / "output.txt"
        f.write_text("data")
        assert fm.skip_if_exists(f, force=False) is True

    def test_skip_if_exists_returns_false_when_file_missing(self, tmp_path):
        fm = FileManager(make_config(tmp_path))
        f = tmp_path / "missing.txt"
        assert fm.skip_if_exists(f, force=False) is False

    def test_skip_if_exists_returns_false_when_force_true(self, tmp_path):
        fm = FileManager(make_config(tmp_path))
        f = tmp_path / "output.txt"
        f.write_text("data")
        assert fm.skip_if_exists(f, force=True) is False

    def test_skip_if_exists_returns_false_when_missing_and_force(self, tmp_path):
        fm = FileManager(make_config(tmp_path))
        f = tmp_path / "missing.txt"
        assert fm.skip_if_exists(f, force=True) is False
