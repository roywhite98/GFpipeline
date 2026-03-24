"""Top-level CLI entry point for gfpipeline."""

from __future__ import annotations

import logging
import sys

import click

from gfpipeline.core.exceptions import ConfigError, PipelineError
from gfpipeline.config.loader import load_config


def setup_logging(verbose: bool) -> None:
    """Configure root logger.

    Args:
        verbose: If True, set DEBUG level; otherwise INFO.
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )


def load_config_from_ctx(ctx: click.Context):
    """Load PipelineConfig from the path stored in ctx.obj['config'].

    Args:
        ctx: Click context carrying global options.

    Returns:
        PipelineConfig instance.

    Raises:
        ConfigError: If the config file does not exist or is invalid.
    """
    config_path = ctx.obj["config"]
    try:
        return load_config(config_path)
    except ConfigError:
        raise
    except FileNotFoundError:
        raise ConfigError(f"Configuration file not found: {config_path}")


@click.group()
@click.option("--config", "-c", default="config.yaml", help="Config file path")
@click.option("--dry-run", is_flag=True, default=False, help="Print commands without executing")
@click.option("--verbose", is_flag=True, default=False, help="Verbose logging")
@click.option("--force", is_flag=True, default=False, help="Force regenerate intermediate files")
@click.pass_context
def cli(ctx: click.Context, config: str, dry_run: bool, verbose: bool, force: bool) -> None:
    """gene-family-pipeline: Gene family analysis pipeline."""
    ctx.ensure_object(dict)
    ctx.obj["config"] = config
    ctx.obj["dry_run"] = dry_run
    ctx.obj["verbose"] = verbose
    ctx.obj["force"] = force
    setup_logging(verbose)


from gfpipeline.cli import cmd_genome_db  # noqa: F401 — registers genome-db subcommand
from gfpipeline.cli import cmd_run  # noqa: F401 — registers run subcommand
from gfpipeline.cli import cmd_stages  # noqa: F401 — registers stage subcommands
from gfpipeline.cli import cmd_refine  # noqa: F401 — registers refine subcommands


def main() -> None:
    """Entry point wrapper that catches PipelineError and exits with code 1."""
    try:
        cli(standalone_mode=False)
    except PipelineError as exc:
        click.echo(f"[ERROR] {exc}", err=True)
        sys.exit(1)
    except SystemExit:
        raise
    except Exception as exc:  # noqa: BLE001
        click.echo(f"[ERROR] Unexpected error: {exc}", err=True)
        sys.exit(1)
