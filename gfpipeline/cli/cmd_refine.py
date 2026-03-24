"""CLI subcommands for the interactive refinement stage.

Provides: refine, refine-tree, refine-domain, refine-motif
"""

from __future__ import annotations

import click

from gfpipeline.cli.main import cli, load_config_from_ctx
from gfpipeline.core.exceptions import ConfigError
from gfpipeline.core.file_manager import FileManager
from gfpipeline.core.runner import ToolRunner
from gfpipeline.stages.refine import RefineStage


def _make_refine_stage(ctx: click.Context) -> tuple[RefineStage, "PipelineConfig"]:  # noqa: F821
    """Load config and instantiate RefineStage from ctx.obj."""
    config = load_config_from_ctx(ctx)
    if config.refinement is None or not config.refinement.idlist:
        raise ConfigError(
            "Missing required config field 'refinement.idlist'. "
            "Please add a 'refinement' section with 'idlist' to your config file."
        )
    runner = ToolRunner(dry_run=ctx.obj["dry_run"], verbose=ctx.obj["verbose"])
    fm = FileManager(config)
    fm.ensure_dirs()
    return RefineStage(config, runner, fm), config


@cli.command("refine")
@click.pass_context
def refine_cmd(ctx: click.Context) -> None:
    """Run the full refine pipeline: extract sequences → tree → domain → motif."""
    stage, _ = _make_refine_stage(ctx)
    stage.run(force=ctx.obj["force"])


@cli.command("refine-tree")
@click.pass_context
def refine_tree_cmd(ctx: click.Context) -> None:
    """Run the refine-tree sub-stage only."""
    stage, _ = _make_refine_stage(ctx)
    stage.run_tree(force=ctx.obj["force"])


@cli.command("refine-domain")
@click.pass_context
def refine_domain_cmd(ctx: click.Context) -> None:
    """Run the refine-domain sub-stage only."""
    stage, _ = _make_refine_stage(ctx)
    stage.run_domain(force=ctx.obj["force"])


@cli.command("refine-motif")
@click.pass_context
def refine_motif_cmd(ctx: click.Context) -> None:
    """Run the refine-motif sub-stage only."""
    stage, _ = _make_refine_stage(ctx)
    stage.run_motif(force=ctx.obj["force"])
