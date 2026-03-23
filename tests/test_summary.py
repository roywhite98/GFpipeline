"""Tests for Summary_Reporter."""

from __future__ import annotations

from pathlib import Path

import pytest

from gfpipeline.config.schema import DatabasesConfig, PipelineConfig
from gfpipeline.core.file_manager import FileManager
from gfpipeline.report.summary import Summary_Reporter


# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------

def _make_config(tmp_path: Path) -> PipelineConfig:
    return PipelineConfig(
        project_name="TESTPROJ",
        data_dir=str(tmp_path / "data"),
        result_dir=str(tmp_path / "results"),
        databases=DatabasesConfig(
            ref_fa=str(tmp_path / "ref.fa"),
            genome=str(tmp_path / "genome.fa"),
            gff3=str(tmp_path / "anno.gff3"),
        ),
    )


def _make_reporter(tmp_path: Path) -> Summary_Reporter:
    config = _make_config(tmp_path)
    fm = FileManager(config)
    fm.ensure_dirs()
    return Summary_Reporter(config, fm)


# ---------------------------------------------------------------------------
# 1. write creates summary.txt
# ---------------------------------------------------------------------------

def test_write_creates_summary_file(tmp_path):
    reporter = _make_reporter(tmp_path)
    out = reporter.write()
    assert out.exists()
    assert out.name == "TESTPROJ.summary.txt"


# ---------------------------------------------------------------------------
# 2. Summary contains project name
# ---------------------------------------------------------------------------

def test_write_contains_project_name(tmp_path):
    reporter = _make_reporter(tmp_path)
    out = reporter.write()
    content = out.read_text(encoding="utf-8")
    assert "TESTPROJ" in content


# ---------------------------------------------------------------------------
# 3. Summary marks missing stage files as "未执行"
# ---------------------------------------------------------------------------

def test_write_marks_missing_files_as_not_executed(tmp_path):
    reporter = _make_reporter(tmp_path)
    out = reporter.write()
    content = out.read_text(encoding="utf-8")
    assert "未执行" in content


# ---------------------------------------------------------------------------
# 4. Summary shows file path when file exists
# ---------------------------------------------------------------------------

def test_write_shows_path_when_file_exists(tmp_path):
    reporter = _make_reporter(tmp_path)
    # Create the identify candidates pep file
    pep_path = reporter.fm.result("identify", "candidates.pep", "fa")
    pep_path.write_text(">Gene1\nMACK\n", encoding="utf-8")

    out = reporter.write()
    content = out.read_text(encoding="utf-8")
    assert str(pep_path) in content


# ---------------------------------------------------------------------------
# 5. Summary shows gene count for idlist files
# ---------------------------------------------------------------------------

def test_write_shows_gene_count_for_idlist(tmp_path):
    reporter = _make_reporter(tmp_path)
    idlist_path = reporter.fm.result("identify", "candidates.gene", "idlist")
    idlist_path.write_text("Gene1\nGene2\nGene3\n", encoding="utf-8")

    out = reporter.write()
    content = out.read_text(encoding="utf-8")
    assert "3 个基因" in content


# ---------------------------------------------------------------------------
# 6. collect returns dict with all expected stage keys
# ---------------------------------------------------------------------------

def test_collect_returns_all_stage_keys(tmp_path):
    reporter = _make_reporter(tmp_path)
    data = reporter.collect()
    expected_stages = {
        "identify", "tree", "domain", "domain-filter",
        "motif", "collinearity", "properties", "trans",
    }
    assert set(data.keys()) == expected_stages


# ---------------------------------------------------------------------------
# 7. Summary contains all stage section headers
# ---------------------------------------------------------------------------

def test_write_contains_all_stage_headers(tmp_path):
    reporter = _make_reporter(tmp_path)
    out = reporter.write()
    content = out.read_text(encoding="utf-8")

    expected_headers = [
        "--- identify 阶段 ---",
        "--- tree 阶段 ---",
        "--- domain 阶段 ---",
        "--- domain-filter 阶段 ---",
        "--- motif 阶段 ---",
        "--- collinearity 阶段 ---",
        "--- properties 阶段 ---",
        "--- trans 阶段 ---",
    ]
    for header in expected_headers:
        assert header in content, f"Missing header: {header}"


# ---------------------------------------------------------------------------
# 8. write returns the correct path
# ---------------------------------------------------------------------------

def test_write_returns_correct_path(tmp_path):
    reporter = _make_reporter(tmp_path)
    out = reporter.write()
    expected = reporter.fm.result_dir / "TESTPROJ.summary.txt"
    assert out == expected


# ---------------------------------------------------------------------------
# Bonus: idlist with empty lines counts only non-empty
# ---------------------------------------------------------------------------

def test_idlist_count_ignores_empty_lines(tmp_path):
    reporter = _make_reporter(tmp_path)
    idlist_path = reporter.fm.result("domain-filter", "candidates", "idlist")
    idlist_path.write_text("Gene1\n\nGene2\n\n", encoding="utf-8")

    out = reporter.write()
    content = out.read_text(encoding="utf-8")
    assert "2 个基因" in content


# ---------------------------------------------------------------------------
# Bonus: TSV row count for collinearity blocks
# ---------------------------------------------------------------------------

def test_tsv_row_count_for_blocks(tmp_path):
    reporter = _make_reporter(tmp_path)
    blocks_path = reporter.fm.result("collinearity", "blocks", "tsv")
    blocks_path.write_text("col1\tcol2\nA\t1\nB\t2\nC\t3\n", encoding="utf-8")

    out = reporter.write()
    content = out.read_text(encoding="utf-8")
    assert "3 个块" in content
