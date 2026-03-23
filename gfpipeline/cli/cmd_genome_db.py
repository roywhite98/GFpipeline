"""CLI subcommand group: genome-db.

Provides blast, gff-qc, gene-index, rep-index, and query sub-subcommands.
When invoked without a sub-subcommand, runs all four build steps in order.
"""

from __future__ import annotations

from pathlib import Path

import click

from gfpipeline.cli.main import cli, load_config_from_ctx
from gfpipeline.core.runner import ToolRunner
from gfpipeline.core.file_manager import FileManager
from gfpipeline.genome_db.blast_db import BlastDbBuilder
from gfpipeline.genome_db.gff_qc import GffQc
from gfpipeline.genome_db.gene_index import GeneIndexBuilder
from gfpipeline.genome_db.rep_index import RepIndexBuilder


@cli.group("genome-db", invoke_without_command=True)
@click.pass_context
def genome_db(ctx: click.Context) -> None:
    """Genome database preparation (blast, gff-qc, gene-index, rep-index)."""
    if ctx.invoked_subcommand is None:
        ctx.invoke(blast_cmd)
        ctx.invoke(gff_qc_cmd)
        ctx.invoke(gene_index_cmd)
        ctx.invoke(rep_index_cmd)


@genome_db.command("blast")
@click.pass_context
def blast_cmd(ctx: click.Context) -> None:
    """Build BLAST databases."""
    config = load_config_from_ctx(ctx)
    force = ctx.obj["force"]
    dry_run = ctx.obj["dry_run"]
    verbose = ctx.obj["verbose"]
    runner = ToolRunner(dry_run=dry_run, verbose=verbose)
    fm = FileManager(config)
    BlastDbBuilder(config, runner, fm).run(force=force)


@genome_db.command("gff-qc")
@click.pass_context
def gff_qc_cmd(ctx: click.Context) -> None:
    """Run GFF3 quality control and correction."""
    config = load_config_from_ctx(ctx)
    force = ctx.obj["force"]
    dry_run = ctx.obj["dry_run"]
    verbose = ctx.obj["verbose"]
    runner = ToolRunner(dry_run=dry_run, verbose=verbose)
    GffQc(config, runner).run(force=force)


@genome_db.command("gene-index")
@click.pass_context
def gene_index_cmd(ctx: click.Context) -> None:
    """Build gene/transcript sequence indices."""
    config = load_config_from_ctx(ctx)
    force = ctx.obj["force"]
    dry_run = ctx.obj["dry_run"]
    verbose = ctx.obj["verbose"]
    runner = ToolRunner(dry_run=dry_run, verbose=verbose)
    fm = FileManager(config)
    GeneIndexBuilder(config, runner, fm).run(force=force)


@genome_db.command("rep-index")
@click.pass_context
def rep_index_cmd(ctx: click.Context) -> None:
    """Build representative transcript index."""
    config = load_config_from_ctx(ctx)
    force = ctx.obj["force"]
    fm = FileManager(config)
    RepIndexBuilder(config, fm).run(force=force)


@genome_db.command("query")
@click.option("--id", "id_", required=True, help="Gene or transcript ID")
@click.option(
    "--type", "seq_type",
    type=click.Choice(["cds", "pep", "genome"]),
    required=True,
    help="Sequence type to retrieve",
)
@click.option("--output", "-o", default=None, help="Output file (default: stdout)")
@click.pass_context
def query_cmd(ctx: click.Context, id_: str, seq_type: str, output: str | None) -> None:
    """Query a sequence by ID."""
    config = load_config_from_ctx(ctx)
    dry_run = ctx.obj["dry_run"]
    verbose = ctx.obj["verbose"]
    runner = ToolRunner(dry_run=dry_run, verbose=verbose)
    fm = FileManager(config)
    result = GeneIndexBuilder(config, runner, fm).query(id_, seq_type)
    if output:
        Path(output).write_text(result)
    else:
        click.echo(result, nl=False)
