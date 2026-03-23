"""CLI subcommand: run.

Runs the full pipeline in order:
    identify → tree → domain → domain-filter → motif → collinearity
Auto-builds the BLAST database if it doesn't exist, then generates a summary report.
"""

from __future__ import annotations

from pathlib import Path

import click

from gfpipeline.cli.main import cli, load_config_from_ctx
from gfpipeline.core.runner import ToolRunner
from gfpipeline.core.file_manager import FileManager
from gfpipeline.stages.identify import IdentifyStage
from gfpipeline.stages.tree import TreeStage
from gfpipeline.stages.domain import DomainStage
from gfpipeline.stages.domain_filter import DomainFilterStage
from gfpipeline.stages.motif import MotifStage
from gfpipeline.stages.collinearity import CollinearityStage
from gfpipeline.genome_db.blast_db import BlastDbBuilder
from gfpipeline.report.summary import Summary_Reporter


def _blast_db_exists(prefix: str) -> bool:
    """Check if protein BLAST DB files exist."""
    return any(Path(prefix + ext).exists() for ext in [".phr", ".pin", ".psq"])


@cli.command("run")
@click.pass_context
def run_cmd(ctx: click.Context) -> None:
    """Run the full pipeline: identify → tree → domain → domain-filter → motif → collinearity."""
    config = load_config_from_ctx(ctx)
    force = ctx.obj["force"]
    dry_run = ctx.obj["dry_run"]
    verbose = ctx.obj["verbose"]

    runner = ToolRunner(dry_run=dry_run, verbose=verbose)
    fm = FileManager(config)
    fm.ensure_dirs()

    # Auto-run genome-db blast if BLAST DB doesn't exist
    if not _blast_db_exists(fm.blast_db_prefix):
        click.echo("BLAST database not found; running genome-db blast first...")
        BlastDbBuilder(config, runner, fm).run(force=force)

    # Run stages in order
    IdentifyStage(config, runner, fm).run(force=force)
    TreeStage(config, runner, fm).run(force=force)
    DomainStage(config, fm).run(force=force)
    DomainFilterStage(config, fm).run(force=force)
    MotifStage(config, runner, fm).run(force=force)
    CollinearityStage(config, runner, fm).run(force=force)

    # Generate summary report
    Summary_Reporter(config, fm).write()
    click.echo("Pipeline complete.")
