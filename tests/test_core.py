"""Unit tests for core layer: exceptions, logger, runner."""

import subprocess
import sys
import logging
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest

from gfpipeline.core.exceptions import (
    PipelineError,
    ConfigError,
    ToolNotFoundError,
    StageInputError,
    ExternalToolError,
    ApiError,
)
from gfpipeline.core.runner import ToolRunner
from gfpipeline.core.logger import setup_logging, get_logger


# ---------------------------------------------------------------------------
# Exception hierarchy tests
# ---------------------------------------------------------------------------

class TestExceptionHierarchy:
    def test_config_error_is_pipeline_error(self):
        assert issubclass(ConfigError, PipelineError)

    def test_tool_not_found_is_pipeline_error(self):
        assert issubclass(ToolNotFoundError, PipelineError)

    def test_stage_input_error_is_pipeline_error(self):
        assert issubclass(StageInputError, PipelineError)

    def test_external_tool_error_is_pipeline_error(self):
        assert issubclass(ExternalToolError, PipelineError)

    def test_api_error_is_pipeline_error(self):
        assert issubclass(ApiError, PipelineError)

    def test_external_tool_error_message_contains_tool_name(self):
        err = ExternalToolError(tool="blastp", cmd=["blastp", "-h"], returncode=1, stderr="bad")
        assert "blastp" in str(err)
        assert "1" in str(err)

    def test_api_error_message_contains_status_code(self):
        err = ApiError(status_code=503, reason="Service Unavailable")
        assert "503" in str(err)

    def test_api_error_without_status_code(self):
        err = ApiError(reason="timeout")
        assert "timeout" in str(err)


# ---------------------------------------------------------------------------
# ToolRunner tests
# ---------------------------------------------------------------------------

class TestToolRunner:
    def test_dry_run_does_not_execute(self):
        runner = ToolRunner(dry_run=True)
        result = runner.run(["false"])  # 'false' always exits 1
        assert result.returncode == 0

    def test_dry_run_shell_does_not_execute(self):
        runner = ToolRunner(dry_run=True)
        result = runner.run_shell("exit 1")
        assert result.returncode == 0

    def test_successful_command_returns_completed_process(self):
        runner = ToolRunner()
        result = runner.run(["echo", "hello"])
        assert result.returncode == 0

    def test_failed_command_raises_external_tool_error(self):
        runner = ToolRunner()
        with pytest.raises(ExternalToolError) as exc_info:
            runner.run(["false"])
        assert exc_info.value.returncode != 0

    def test_failed_shell_command_raises_external_tool_error(self):
        runner = ToolRunner()
        with pytest.raises(ExternalToolError):
            runner.run_shell("exit 42")

    def test_verbose_runner_logs_command(self, caplog):
        runner = ToolRunner(verbose=True)
        with caplog.at_level(logging.DEBUG, logger="gfpipeline.core.runner"):
            runner.run(["echo", "test"])
        assert any("echo" in r.message for r in caplog.records)

    def test_run_with_cwd(self, tmp_path):
        runner = ToolRunner()
        result = runner.run(["pwd"], cwd=str(tmp_path))
        assert result.returncode == 0

    def test_external_tool_error_stores_attributes(self):
        runner = ToolRunner()
        with pytest.raises(ExternalToolError) as exc_info:
            runner.run(["false"])
        err = exc_info.value
        assert err.tool == "false"
        assert err.returncode != 0


# ---------------------------------------------------------------------------
# Logger tests
# ---------------------------------------------------------------------------

class TestLogger:
    def test_setup_logging_creates_log_file(self, tmp_path):
        setup_logging(result_dir=tmp_path, verbose=False)
        log = get_logger("test.setup")
        log.info("test message")
        log_file = tmp_path / "pipeline.log"
        assert log_file.exists()
        assert "test message" in log_file.read_text()

    def test_setup_logging_verbose_sets_debug_level(self, tmp_path):
        setup_logging(result_dir=tmp_path, verbose=True)
        root = logging.getLogger()
        assert root.level == logging.DEBUG

    def test_setup_logging_default_sets_info_level(self, tmp_path):
        setup_logging(result_dir=tmp_path, verbose=False)
        root = logging.getLogger()
        assert root.level == logging.INFO

    def test_get_logger_returns_named_logger(self, tmp_path):
        setup_logging(result_dir=tmp_path)
        logger = get_logger("gfpipeline.test")
        assert logger.name == "gfpipeline.test"

    def test_log_file_appends_across_calls(self, tmp_path):
        setup_logging(result_dir=tmp_path)
        get_logger("t").info("first")
        setup_logging(result_dir=tmp_path)
        get_logger("t").info("second")
        content = (tmp_path / "pipeline.log").read_text()
        assert "first" in content
        assert "second" in content

    def test_setup_logging_creates_result_dir(self, tmp_path):
        nested = tmp_path / "nested" / "results"
        setup_logging(result_dir=nested)
        assert (nested / "pipeline.log").exists()
