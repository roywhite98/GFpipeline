"""Tests for DomainStage."""

from __future__ import annotations

import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from gfpipeline.config.schema import DatabasesConfig, DomainConfig, PipelineConfig
from gfpipeline.core.exceptions import ApiError, StageInputError
from gfpipeline.core.file_manager import FileManager
from gfpipeline.stages.domain import DomainStage


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_config(tmp_path: Path, domain_cfg: DomainConfig | None = None) -> PipelineConfig:
    return PipelineConfig(
        project_name="TEST",
        data_dir=str(tmp_path / "data"),
        result_dir=str(tmp_path / "results"),
        databases=DatabasesConfig(
            ref_fa=str(tmp_path / "ref.fa"),
            genome=str(tmp_path / "genome.fa"),
            gff3=str(tmp_path / "anno.gff3"),
        ),
        domain=domain_cfg or DomainConfig(),
    )


def _make_stage(tmp_path: Path) -> DomainStage:
    config = _make_config(tmp_path)
    fm = FileManager(config)
    return DomainStage(config, fm)


def _fake_submit_response() -> str:
    return textwrap.dedent("""\
        #Batch CD-search tool
        #cdsid    QM3-qcdsearch-12345678-90AB
        #datatype    hitsFull Results
        #status    3
    """)


def _fake_status_inprogress() -> str:
    return textwrap.dedent("""\
        #status    3
        #cdsid    QM3-qcdsearch-12345678-90AB
    """)


def _fake_status_complete() -> str:
    return textwrap.dedent("""\
        #status    0
        #cdsid    QM3-qcdsearch-12345678-90AB
    """)


def _fake_result_text() -> str:
    return textwrap.dedent("""\
        #Batch CD-search tool
        #cdsid    QM3-qcdsearch-12345678-90AB
        Query\tHit_type\tPSSM_ID\tFrom\tTo\tE-Value\tBitscore\tAccession\tShort_name\tIncomplete\tSuperfamily
        gene1\tspecific\t12345\t1\t100\t1e-10\t200\tcd12345\tDomain_A\t-\tcl99999
    """)


# ---------------------------------------------------------------------------
# Tests: _parse_rid
# ---------------------------------------------------------------------------

def test_parse_rid_success(tmp_path):
    stage = _make_stage(tmp_path)
    rid = stage._parse_rid(_fake_submit_response())
    assert rid == "QM3-qcdsearch-12345678-90AB"


def test_parse_rid_missing_raises(tmp_path):
    stage = _make_stage(tmp_path)
    with pytest.raises(ApiError):
        stage._parse_rid("no cdsid here")


# ---------------------------------------------------------------------------
# Tests: _parse_status
# ---------------------------------------------------------------------------

def test_parse_status_complete(tmp_path):
    stage = _make_stage(tmp_path)
    assert stage._parse_status(_fake_status_complete()) == 0


def test_parse_status_inprogress(tmp_path):
    stage = _make_stage(tmp_path)
    assert stage._parse_status(_fake_status_inprogress()) == 3


def test_parse_status_missing_returns_3(tmp_path):
    stage = _make_stage(tmp_path)
    assert stage._parse_status("no status line") == 3


# ---------------------------------------------------------------------------
# Tests: submit_batch
# ---------------------------------------------------------------------------

def test_submit_batch_returns_rid(tmp_path):
    stage = _make_stage(tmp_path)
    pep_fa = tmp_path / "test.fa"
    pep_fa.write_text(">gene1\nMACDEFGH\n")

    mock_resp = MagicMock()
    mock_resp.status_code = 200
    mock_resp.text = _fake_submit_response()

    with patch("gfpipeline.stages.domain.requests.post", return_value=mock_resp):
        rid = stage.submit_batch(pep_fa)

    assert rid == "QM3-qcdsearch-12345678-90AB"


def test_submit_batch_non200_raises(tmp_path):
    stage = _make_stage(tmp_path)
    pep_fa = tmp_path / "test.fa"
    pep_fa.write_text(">gene1\nMACDEFGH\n")

    mock_resp = MagicMock()
    mock_resp.status_code = 503
    mock_resp.reason = "Service Unavailable"

    with patch("gfpipeline.stages.domain.requests.post", return_value=mock_resp):
        with pytest.raises(ApiError) as exc_info:
            stage.submit_batch(pep_fa)
    assert exc_info.value.status_code == 503


# ---------------------------------------------------------------------------
# Tests: poll_status
# ---------------------------------------------------------------------------

def test_poll_status_completes_immediately(tmp_path):
    stage = _make_stage(tmp_path)

    mock_resp = MagicMock()
    mock_resp.status_code = 200
    mock_resp.text = _fake_status_complete()

    with patch("gfpipeline.stages.domain.requests.post", return_value=mock_resp):
        # Should return without raising
        stage.poll_status("QM3-qcdsearch-12345678-90AB")


def test_poll_status_error_code_raises(tmp_path):
    stage = _make_stage(tmp_path)

    mock_resp = MagicMock()
    mock_resp.status_code = 200
    mock_resp.text = "#status    2\n#cdsid    QM3-qcdsearch-12345678-90AB\n"

    with patch("gfpipeline.stages.domain.requests.post", return_value=mock_resp):
        with pytest.raises(ApiError):
            stage.poll_status("QM3-qcdsearch-12345678-90AB")


def test_poll_status_non200_raises(tmp_path):
    stage = _make_stage(tmp_path)

    mock_resp = MagicMock()
    mock_resp.status_code = 500
    mock_resp.reason = "Internal Server Error"

    with patch("gfpipeline.stages.domain.requests.post", return_value=mock_resp):
        with pytest.raises(ApiError) as exc_info:
            stage.poll_status("some-rid")
    assert exc_info.value.status_code == 500


def test_poll_status_polls_until_complete(tmp_path):
    """Verify that poll_status keeps polling when status=3 and stops at status=0."""
    stage = _make_stage(tmp_path)

    responses = [
        MagicMock(status_code=200, text="#status    3\n"),
        MagicMock(status_code=200, text="#status    3\n"),
        MagicMock(status_code=200, text="#status    0\n"),
    ]

    with patch("gfpipeline.stages.domain.requests.post", side_effect=responses):
        with patch("gfpipeline.stages.domain.time.sleep"):
            stage.poll_status("some-rid")


# ---------------------------------------------------------------------------
# Tests: download_result
# ---------------------------------------------------------------------------

def test_download_result_saves_file(tmp_path):
    stage = _make_stage(tmp_path)
    output = tmp_path / "result.txt"

    mock_resp = MagicMock()
    mock_resp.status_code = 200
    mock_resp.text = _fake_result_text()

    with patch("gfpipeline.stages.domain.requests.post", return_value=mock_resp):
        stage.download_result("some-rid", output)

    assert output.exists()
    assert "gene1" in output.read_text()


def test_download_result_non200_raises(tmp_path):
    stage = _make_stage(tmp_path)
    output = tmp_path / "result.txt"

    mock_resp = MagicMock()
    mock_resp.status_code = 404
    mock_resp.reason = "Not Found"

    with patch("gfpipeline.stages.domain.requests.post", return_value=mock_resp):
        with pytest.raises(ApiError):
            stage.download_result("some-rid", output)


# ---------------------------------------------------------------------------
# Tests: run
# ---------------------------------------------------------------------------

def test_run_raises_if_pep_missing(tmp_path):
    stage = _make_stage(tmp_path)
    with pytest.raises(StageInputError):
        stage.run()


def test_run_skips_if_output_exists(tmp_path):
    config = _make_config(tmp_path)
    fm = FileManager(config)
    stage = DomainStage(config, fm)

    # Create required input
    fm.ensure_dirs()
    stage._pep_fa.write_text(">gene1\nMACDEFGH\n")
    # Create output so it gets skipped
    stage._cdd_output.write_text("existing content")

    with patch("gfpipeline.stages.domain.requests.post") as mock_post:
        stage.run(force=False)
        mock_post.assert_not_called()


def test_run_full_flow(tmp_path):
    """Integration test: run() submits, polls, downloads, and writes output."""
    config = _make_config(tmp_path)
    fm = FileManager(config)
    fm.ensure_dirs()
    stage = DomainStage(config, fm)

    # Create input pep.fa
    stage._pep_fa.write_text(">gene1\nMACDEFGH\n>gene2\nMKLMNOP\n")

    submit_resp = MagicMock(status_code=200, text=_fake_submit_response())
    poll_resp = MagicMock(status_code=200, text=_fake_status_complete())
    download_resp = MagicMock(status_code=200, text=_fake_result_text())

    with patch(
        "gfpipeline.stages.domain.requests.post",
        side_effect=[submit_resp, poll_resp, download_resp],
    ):
        stage.run()

    assert stage._cdd_output.exists()
    content = stage._cdd_output.read_text()
    assert "gene1" in content
