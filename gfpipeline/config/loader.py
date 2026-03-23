"""YAML configuration loader with Pydantic validation."""

from __future__ import annotations

from pathlib import Path
from typing import Union

import yaml
from pydantic import ValidationError

from gfpipeline.core.exceptions import ConfigError
from gfpipeline.config.schema import PipelineConfig


def load_config(path: Union[str, Path]) -> PipelineConfig:
    """Load and validate a YAML configuration file.

    Args:
        path: Path to the YAML configuration file.

    Returns:
        Validated PipelineConfig instance.

    Raises:
        ConfigError: If the file is not found, cannot be parsed, or fails
                     Pydantic validation (missing/invalid fields).
    """
    config_path = Path(path)

    if not config_path.exists():
        raise ConfigError(f"Configuration file not found: {config_path}")

    try:
        with config_path.open("r", encoding="utf-8") as fh:
            raw = yaml.safe_load(fh)
    except yaml.YAMLError as exc:
        raise ConfigError(f"Failed to parse YAML file '{config_path}': {exc}") from exc

    if not isinstance(raw, dict):
        raise ConfigError(
            f"Configuration file '{config_path}' must contain a YAML mapping, "
            f"got {type(raw).__name__}"
        )

    try:
        return PipelineConfig(**raw)
    except ValidationError as exc:
        missing = []
        for error in exc.errors():
            loc = ".".join(str(part) for part in error["loc"])
            missing.append(loc)
        fields_str = ", ".join(missing)
        raise ConfigError(
            f"Configuration validation failed. Missing or invalid fields: {fields_str}"
        ) from exc
