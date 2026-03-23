"""Logging module for gfpipeline.

Writes to both stderr and Result_Dir/pipeline.log.
Supports INFO/DEBUG level switching via setup_logging().
"""

import logging
import sys
from pathlib import Path

_FORMATTER = logging.Formatter(
    fmt="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

_file_handler: logging.FileHandler | None = None


def setup_logging(result_dir: str | Path, verbose: bool = False) -> None:
    """Configure root logger to write to stderr and Result_Dir/pipeline.log.

    Args:
        result_dir: Path to the result directory where pipeline.log will be written.
        verbose: If True, set level to DEBUG; otherwise INFO.
    """
    global _file_handler

    level = logging.DEBUG if verbose else logging.INFO
    root = logging.getLogger()
    root.setLevel(level)

    # Remove existing handlers to avoid duplicate output on repeated calls
    root.handlers.clear()

    # stderr handler
    stderr_handler = logging.StreamHandler(sys.stderr)
    stderr_handler.setLevel(level)
    stderr_handler.setFormatter(_FORMATTER)
    root.addHandler(stderr_handler)

    # File handler
    log_path = Path(result_dir) / "pipeline.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    _file_handler = logging.FileHandler(log_path, mode="a", encoding="utf-8")
    _file_handler.setLevel(level)
    _file_handler.setFormatter(_FORMATTER)
    root.addHandler(_file_handler)


def get_logger(name: str) -> logging.Logger:
    """Return a named logger.

    Args:
        name: Logger name, typically __name__ of the calling module.
    """
    return logging.getLogger(name)
