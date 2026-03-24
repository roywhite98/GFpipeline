"""CLI subcommands for individual pipeline stages."""

from __future__ import annotations

import click

from gfpipeline.cli.main import cli, load_config_from_ctx
from gfpipeline.core.runner import ToolRunner
from gfpipeline.core.file_manager import FileManager
from gfpipeline.core.tool_checker import check_tools
from gfpipeline.stages.identify import IdentifyStage
from gfpipeline.stages.tree import TreeStage
from gfpipeline.stages.domain import DomainStage
from gfpipeline.stages.domain_filter import DomainFilterStage
from gfpipeline.stages.motif import MotifStage
from gfpipeline.stages.collinearity import CollinearityStage
from gfpipeline.stages.properties import PropertiesStage
from gfpipeline.stages.trans import TransStage


@cli.command("identify")
@click.pass_context
def identify_cmd(ctx: click.Context) -> None:
    """Run the identify stage."""
    config = load_config_from_ctx(ctx)
    force = ctx.obj["force"]
    dry_run = ctx.obj["dry_run"]
    verbose = ctx.obj["verbose"]
    check_tools("identify", config)
    runner = ToolRunner(dry_run=dry_run, verbose=verbose)
    fm = FileManager(config)
    IdentifyStage(config, runner, fm).run(force=force)


@cli.command("tree")
@click.pass_context
def tree_cmd(ctx: click.Context) -> None:
    """Run the tree stage."""
    config = load_config_from_ctx(ctx)
    force = ctx.obj["force"]
    dry_run = ctx.obj["dry_run"]
    verbose = ctx.obj["verbose"]
    check_tools("tree", config)
    runner = ToolRunner(dry_run=dry_run, verbose=verbose)
    fm = FileManager(config)
    TreeStage(config, runner, fm).run(force=force)


@cli.command("domain")
@click.pass_context
def domain_cmd(ctx: click.Context) -> None:
    """Run the domain stage."""
    config = load_config_from_ctx(ctx)
    force = ctx.obj["force"]
    fm = FileManager(config)
    DomainStage(config, fm).run(force=force)


@cli.command("domain-filter")
@click.pass_context
def domain_filter_cmd(ctx: click.Context) -> None:
    """Run the domain-filter stage."""
    config = load_config_from_ctx(ctx)
    force = ctx.obj["force"]
    fm = FileManager(config)
    DomainFilterStage(config, fm).run(force=force)


@cli.command("motif")
@click.pass_context
def motif_cmd(ctx: click.Context) -> None:
    """Run the motif stage."""
    config = load_config_from_ctx(ctx)
    force = ctx.obj["force"]
    dry_run = ctx.obj["dry_run"]
    verbose = ctx.obj["verbose"]
    check_tools("motif", config)
    runner = ToolRunner(dry_run=dry_run, verbose=verbose)
    fm = FileManager(config)
    MotifStage(config, runner, fm).run(force=force)


@cli.command("motif-filter")
@click.pass_context
def motif_filter_cmd(ctx: click.Context) -> None:
    """Re-run motif filtering only (skip MEME/FIMO, use existing fimo.tsv)."""
    config = load_config_from_ctx(ctx)
    fm = FileManager(config)
    runner = ToolRunner(dry_run=ctx.obj["dry_run"], verbose=ctx.obj["verbose"])
    MotifStage(config, runner, fm).run_filter_only()


@cli.command("collinearity")
@click.pass_context
def collinearity_cmd(ctx: click.Context) -> None:
    """Run the collinearity stage."""
    config = load_config_from_ctx(ctx)
    force = ctx.obj["force"]
    dry_run = ctx.obj["dry_run"]
    verbose = ctx.obj["verbose"]
    check_tools("collinearity", config)
    runner = ToolRunner(dry_run=dry_run, verbose=verbose)
    fm = FileManager(config)
    CollinearityStage(config, runner, fm).run(force=force)


@cli.command("properties")
@click.pass_context
def properties_cmd(ctx: click.Context) -> None:
    """Run the properties stage."""
    config = load_config_from_ctx(ctx)
    force = ctx.obj["force"]
    fm = FileManager(config)
    PropertiesStage(config, fm).run(force=force)


@cli.command("trans")
@click.pass_context
def trans_cmd(ctx: click.Context) -> None:
    """Run the trans stage."""
    config = load_config_from_ctx(ctx)
    force = ctx.obj["force"]
    fm = FileManager(config)
    TransStage(config, fm).run(force=force)
