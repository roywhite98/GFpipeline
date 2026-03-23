"""BLAST database builder for genome-db stage.

Builds protein and/or nucleotide BLAST databases using makeblastdb.
"""

from __future__ import annotations

from pathlib import Path

from gfpipeline.config.schema import PipelineConfig
from gfpipeline.core.exceptions import PipelineError, StageInputError
from gfpipeline.core.file_manager import FileManager
from gfpipeline.core.logger import get_logger
from gfpipeline.core.runner import ToolRunner

log = get_logger(__name__)

# File extensions that indicate a BLAST database exists
_PROT_EXTS = [".phr", ".pin", ".psq"]
_NUCL_EXTS = [".nhr", ".nin", ".nsq"]


def _db_exists(prefix: str, exts: list[str]) -> bool:
    """Return True if all expected BLAST db files exist for the given prefix."""
    return all(Path(prefix + ext).exists() for ext in exts)


class BlastDbBuilder:
    """Build BLAST databases for protein and/or nucleotide sequences.

    Args:
        config: Pipeline configuration.
        runner: External tool runner.
        fm: File manager (unused directly but kept for API consistency).
    """

    def __init__(self, config: PipelineConfig, runner: ToolRunner, fm: FileManager) -> None:
        self.config = config
        self.runner = runner
        self.fm = fm

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(self, force: bool = False) -> None:
        """Build BLAST databases according to configuration.

        Steps:
        1. Validate that at least one build type is enabled.
        2. Check input files exist.
        3. Build protein DB from rep_pep if build_prot=True.
        4. Build nucleotide DB if build_nucl=True.
        5. Log statistics.

        Args:
            force: If True, rebuild even if DB files already exist.

        Raises:
            PipelineError: If both build_prot and build_nucl are False.
            StageInputError: If required input files are missing.
        """
        cfg = self.config.genome_db
        db_cfg = self.config.databases

        if not cfg.build_prot and not cfg.build_nucl:
            raise PipelineError(
                "genome-db blast: both build_prot and build_nucl are False. "
                "Enable at least one build type in the configuration."
            )

        if cfg.build_prot:
            self._build_prot(str(self.fm.rep_pep), self.fm.blast_db_prefix, force)

        if cfg.build_nucl:
            self._build_nucl(db_cfg.genome, force)

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _build_prot(self, pep_path: str, blast_db_prefix: str, force: bool) -> None:
        """Build protein BLAST database."""
        pep = Path(pep_path)
        if not pep.exists():
            raise StageInputError(
                f"genome-db blast: protein FASTA not found: {pep}"
            )

        if not force and _db_exists(blast_db_prefix, _PROT_EXTS):
            log.info("Skipping protein BLAST DB (already exists): %s", blast_db_prefix)
            return

        log.info("Building protein BLAST DB: %s -> %s", pep, blast_db_prefix)
        self.runner.run([
            self.config.tools.makeblastdb,
            "-in", str(pep),
            "-dbtype", "prot",
            "-out", blast_db_prefix,
        ])

        self._log_stats(blast_db_prefix, "prot", _PROT_EXTS)

    def _build_nucl(self, genome_path: str, force: bool) -> None:
        """Build nucleotide BLAST database."""
        genome = Path(genome_path)
        if not genome.exists():
            raise StageInputError(
                f"genome-db blast: genome FASTA not found: {genome}"
            )

        nucl_prefix = str(genome) + ".db"

        if not force and _db_exists(nucl_prefix, _NUCL_EXTS):
            log.info("Skipping nucleotide BLAST DB (already exists): %s", nucl_prefix)
            return

        log.info("Building nucleotide BLAST DB: %s -> %s", genome, nucl_prefix)
        self.runner.run([
            self.config.tools.makeblastdb,
            "-in", str(genome),
            "-dbtype", "nucl",
            "-out", nucl_prefix,
        ])

        self._log_stats(nucl_prefix, "nucl", _NUCL_EXTS)

    def _log_stats(self, prefix: str, db_type: str, exts: list[str]) -> None:
        """Log statistics about the built database."""
        files = [prefix + ext for ext in exts if Path(prefix + ext).exists()]
        log.info(
            "BLAST DB built: type=%s, prefix=%s, files=%s",
            db_type, prefix, files,
        )
