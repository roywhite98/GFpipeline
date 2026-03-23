"""Tool dependency checking for gfpipeline stages."""

from __future__ import annotations

import os
import shutil

from gfpipeline.config.schema import PipelineConfig
from gfpipeline.core.exceptions import ToolNotFoundError

# Base tool requirements per stage
STAGE_TOOLS: dict[str, list[str]] = {
    "identify":             ["muscle", "hmmbuild", "hmmsearch", "hmmemit", "blastp"],
    "tree":                 ["muscle", "trimal", "iqtree2"],
    "domain":               [],
    "motif":                ["meme", "fimo"],
    "collinearity":         ["blastp"],   # mcscanx or jcvi added dynamically
    "genome-db blast":      ["makeblastdb"],
    "genome-db gff-qc":     ["samtools"],
    "genome-db gene-index": ["samtools"],
    "genome-db rep-index":  [],
}


def is_executable(tool_path: str) -> bool:
    """Check whether a tool is executable via PATH or as an absolute/relative path."""
    # If it looks like a path (contains a separator), check directly
    if os.sep in tool_path or "/" in tool_path:
        return os.path.isfile(tool_path) and os.access(tool_path, os.X_OK)
    # Otherwise search PATH
    return shutil.which(tool_path) is not None


def check_tools(stage: str, config: PipelineConfig) -> None:
    """Check that all tools required for the given stage are executable.

    Raises ToolNotFoundError listing the first missing tool with install hint.
    """
    tools_to_check = list(STAGE_TOOLS.get(stage, []))

    # Resolve tool names from config
    tool_name_map: dict[str, str] = {
        "muscle":      config.tools.muscle,
        "hmmbuild":    config.tools.hmmbuild,
        "hmmsearch":   config.tools.hmmsearch,
        "hmmemit":     config.tools.hmmemit,
        "blastp":      config.tools.blastp,
        "trimal":      config.tools.trimal,
        "iqtree2":     config.tools.iqtree,
        "meme":        config.tools.meme,
        "fimo":        config.tools.fimo,
        "makeblastdb": config.tools.makeblastdb,
        "samtools":    config.tools.samtools,
        "minimap2":    config.tools.minimap2,
        "stringtie":   config.tools.stringtie,
        "MCScanX":     config.tools.mcscanx,
    }

    # Dynamically add collinearity tool (jcvi or mcscanx)
    if stage == "collinearity":
        if config.collinearity.tool == "mcscanx":
            tools_to_check.append("MCScanX")
        else:
            # jcvi is a Python package; check via python -m jcvi availability
            tools_to_check.append("jcvi")

    # Dynamically add optional gff-qc tools
    if stage == "genome-db gff-qc":
        if config.gff_qc.transcript_fa:
            tools_to_check.append("minimap2")
        if config.gff_qc.rna_bam:
            tools_to_check.append("stringtie")

    missing: list[str] = []
    for tool_key in tools_to_check:
        resolved = tool_name_map.get(tool_key, tool_key)
        if not is_executable(resolved):
            missing.append(resolved)

    if missing:
        tool_list = ", ".join(missing)
        raise ToolNotFoundError(
            f"Required tool(s) not found for stage '{stage}': {tool_list}. "
            f"Please install them and ensure they are in your PATH or configure "
            f"their paths in the config file."
        )
