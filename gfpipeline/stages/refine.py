"""Refine stage: interactive secondary filtering and re-analysis.

Reads a user-provided ID list, extracts protein sequences from Rep_Pep_DB,
then re-runs tree, domain, and motif analyses on the refined sequence set.
"""

from __future__ import annotations

import logging
from pathlib import Path

from gfpipeline.config.schema import PipelineConfig
from gfpipeline.core.exceptions import StageInputError
from gfpipeline.core.file_manager import FileManager, RefineFileManager
from gfpipeline.core.runner import ToolRunner
from gfpipeline.core.sequence import extract_fasta
from gfpipeline.genome_db.gene_index import strip_id_prefix
from gfpipeline.stages.domain import DomainStage
from gfpipeline.stages.motif import MotifStage
from gfpipeline.stages.tree import TreeStage

log = logging.getLogger(__name__)


class RefineStage:
    """Interactive secondary filtering and re-analysis stage.

    Responsibilities:
        1. Read and parse Refinement_IDList
        2. Convert gene IDs to representative transcript IDs (via gene2transcript.tsv)
        3. Extract protein sequences from Rep_Pep_DB → refine.candidates.pep.fa
        4. Schedule refine-tree / refine-domain / refine-motif
    """

    def __init__(self, config: PipelineConfig, runner: ToolRunner, fm: FileManager) -> None:
        self.config = config
        self.runner = runner
        self.fm = fm
        self._rfm = RefineFileManager(config)

    # ------------------------------------------------------------------
    # Path helpers
    # ------------------------------------------------------------------

    @property
    def _candidates_pep(self) -> Path:
        """Path to refine.candidates.pep.fa (via RefineFileManager)."""
        return self._rfm.result("candidates", "pep", "fa")

    # ------------------------------------------------------------------
    # Core methods
    # ------------------------------------------------------------------

    def read_idlist(self, path: Path) -> list[str]:
        """Read idlist file, ignore blank lines and # comment lines, deduplicate.

        Args:
            path: Path to the idlist file.

        Returns:
            Deduplicated list of valid IDs (order-preserving).

        Raises:
            StageInputError: If the file does not exist or no valid IDs remain.
        """
        if not path.exists():
            raise StageInputError(
                f"IDList file '{path}' does not exist. "
                "Please place the idlist file in the data/ directory."
            )

        seen: set[str] = set()
        ids: list[str] = []
        with open(path) as fh:
            for line in fh:
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                if stripped not in seen:
                    seen.add(stripped)
                    ids.append(stripped)

        if not ids:
            raise StageInputError(
                f"IDList file '{path}' contains no valid IDs after filtering "
                "blank lines and comments. Please add gene or transcript IDs."
            )

        log.info("Read %d unique IDs from %s", len(ids), path)
        return ids

    def resolve_ids(self, raw_ids: list[str]) -> list[str]:
        """Convert gene IDs to representative transcript IDs; pass transcript IDs through.

        Reads gene2transcript.tsv from genome_db index_dir. For each input ID:
        - If it maps to a gene in the TSV (is_representative=True), replace with transcript ID.
        - If it looks like a transcript ID (not found as a gene key), pass through unchanged.
        - IDs with no mapping found are logged as WARNING.

        Args:
            raw_ids: List of gene IDs and/or transcript IDs.

        Returns:
            List of resolved transcript IDs (order-preserving, deduplicated).
        """
        index_dir = Path(self.config.genome_db.index_dir)
        g2t_path = index_dir / "gene2transcript.tsv"

        # Build gene_id -> representative transcript_id mapping
        gene_to_rep: dict[str, str] = {}
        # Also collect all known transcript IDs for pass-through detection
        all_transcript_ids: set[str] = set()

        if g2t_path.exists():
            with open(g2t_path) as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("gene_id"):
                        continue
                    parts = line.split("\t")
                    if len(parts) < 3:
                        continue
                    gene_id = strip_id_prefix(parts[0].strip())
                    transcript_id = strip_id_prefix(parts[1].strip())
                    is_rep = parts[2].strip().lower() in ("true", "1", "yes")
                    all_transcript_ids.add(transcript_id)
                    if is_rep:
                        gene_to_rep[gene_id] = transcript_id
        else:
            log.warning("gene2transcript.tsv not found at %s; all IDs treated as transcript IDs", g2t_path)

        resolved: list[str] = []
        seen: set[str] = set()

        for raw_id in raw_ids:
            # Strip known prefixes for lookup
            clean_id = strip_id_prefix(raw_id)

            if clean_id in gene_to_rep:
                # It's a gene ID — map to representative transcript
                t_id = gene_to_rep[clean_id]
                if t_id not in seen:
                    seen.add(t_id)
                    resolved.append(t_id)
            elif clean_id in all_transcript_ids:
                # It's a known transcript ID — pass through
                if clean_id not in seen:
                    seen.add(clean_id)
                    resolved.append(clean_id)
            else:
                # Unknown ID — treat as transcript ID (pass through) but warn
                log.warning("No mapping found for ID '%s'; treating as transcript ID", raw_id)
                if clean_id not in seen:
                    seen.add(clean_id)
                    resolved.append(clean_id)

        log.info("Resolved %d IDs → %d transcript IDs", len(raw_ids), len(resolved))
        return resolved

    def extract_sequences(self, transcript_ids: list[str]) -> Path:
        """Extract sequences from Rep_Pep_DB to refine.candidates.pep.fa.

        Args:
            transcript_ids: List of transcript IDs to extract.

        Returns:
            Path to the output FASTA file.

        Raises:
            StageInputError: If ALL IDs are not found in the database.
        """
        db_path = self.fm.rep_pep
        output_path = self._candidates_pep

        self._rfm.ensure_dirs()

        count = extract_fasta(db_path, transcript_ids, output_path)

        missing = len(transcript_ids) - count
        if missing > 0:
            log.warning(
                "%d/%d IDs not found in Rep_Pep_DB (%s)",
                missing, len(transcript_ids), db_path,
            )

        if count == 0:
            raise StageInputError(
                f"None of the {len(transcript_ids)} IDs were found in Rep_Pep_DB '{db_path}'. "
                "Please check that the IDs match the database."
            )

        log.info(
            "Extracted %d/%d sequences → %s",
            count, len(transcript_ids), output_path,
        )
        return output_path

    # ------------------------------------------------------------------
    # Sub-stage runners
    # ------------------------------------------------------------------

    def run_tree(self, force: bool = False) -> None:
        """Run refine-tree sub-stage using RefineFileManager."""
        if not self._candidates_pep.exists():
            raise StageInputError(
                f"Refined candidates file '{self._candidates_pep}' does not exist. "
                "Please run 'gfpipeline refine' to extract sequences first."
            )
        log.info("Starting refine-tree sub-stage")
        tree_stage = TreeStage(self.config, self.runner, self._rfm)
        # Override the input path: TreeStage looks for identify.candidates.pep.fa,
        # but we need it to use refine.candidates.pep.fa.
        # We monkey-patch _pep_fa via a subclass approach using the rfm.
        # Since TreeStage._pep_fa = fm.result("identify", "candidates.pep", "fa"),
        # and rfm.result("identify", ...) → refine.identify.candidates.pep.fa (wrong),
        # we need a custom approach: create a wrapper that returns the right path.
        _RefinedTreeStage(self.config, self.runner, self._rfm, self._candidates_pep).run(force)
        log.info("refine-tree complete → %s.*", self._rfm.result("tree", "pep.trimed", "afa"))

    def run_domain(self, force: bool = False) -> None:
        """Run refine-domain sub-stage using RefineFileManager."""
        if not self._candidates_pep.exists():
            raise StageInputError(
                f"Refined candidates file '{self._candidates_pep}' does not exist. "
                "Please run 'gfpipeline refine' to extract sequences first."
            )
        log.info("Starting refine-domain sub-stage")
        _RefinedDomainStage(self.config, self._rfm, self._candidates_pep).run(force)
        log.info("refine-domain complete → %s", self._rfm.result("domain", "cdd", "txt"))

    def run_motif(self, force: bool = False) -> None:
        """Run refine-motif sub-stage using RefineFileManager."""
        if not self._candidates_pep.exists():
            raise StageInputError(
                f"Refined candidates file '{self._candidates_pep}' does not exist. "
                "Please run 'gfpipeline refine' to extract sequences first."
            )
        log.info("Starting refine-motif sub-stage")
        _RefinedMotifStage(self.config, self.runner, self._rfm, self._candidates_pep).run(force)
        log.info(
            "refine-motif complete → %s",
            self._rfm.result("motif-filter", "candidates", "idlist"),
        )

    # ------------------------------------------------------------------
    # Main run
    # ------------------------------------------------------------------

    def run(self, force: bool = False) -> None:
        """Execute the full refine pipeline: extract → tree → domain → motif."""
        # Step 1: read and resolve IDs
        idlist_path = Path(self.config.refinement.idlist)
        raw_ids = self.read_idlist(idlist_path)
        transcript_ids = self.resolve_ids(raw_ids)

        # Step 2: extract sequences
        candidates = self.extract_sequences(transcript_ids)
        log.info("Sequence extraction complete → %s", candidates)

        # Step 3: run sub-stages
        self.run_tree(force)
        self.run_domain(force)
        self.run_motif(force)

        log.info("Refine stage complete.")
        log.info("  candidates:  %s", self._candidates_pep)
        log.info("  tree:        %s.*", self._rfm.result("tree", "pep.trimed", "afa"))
        log.info("  domain:      %s", self._rfm.result("domain", "cdd", "txt"))
        log.info("  motif:       %s", self._rfm.result("motif-filter", "candidates", "idlist"))


# ---------------------------------------------------------------------------
# Private helper subclasses that override the input pep path
# ---------------------------------------------------------------------------

class _RefinedTreeStage(TreeStage):
    """TreeStage variant that uses a custom input pep.fa path."""

    def __init__(
        self,
        config: PipelineConfig,
        runner: ToolRunner,
        fm: FileManager,
        pep_fa: Path,
    ) -> None:
        super().__init__(config, runner, fm)
        self.__pep_fa = pep_fa

    @property
    def _pep_fa(self) -> Path:
        return self.__pep_fa


class _RefinedDomainStage(DomainStage):
    """DomainStage variant that uses a custom input pep.fa path."""

    def __init__(
        self,
        config: PipelineConfig,
        fm: FileManager,
        pep_fa: Path,
    ) -> None:
        super().__init__(config, fm)
        self.__pep_fa = pep_fa

    @property
    def _pep_fa(self) -> Path:
        return self.__pep_fa


class _RefinedMotifStage(MotifStage):
    """MotifStage variant that uses a custom input pep.fa path and refine-prefixed dirs."""

    def __init__(
        self,
        config: PipelineConfig,
        runner: ToolRunner,
        fm: FileManager,
        pep_fa: Path,
    ) -> None:
        super().__init__(config, runner, fm)
        self.__pep_fa = pep_fa

    @property
    def _pep_fa(self) -> Path:
        return self.__pep_fa

    @property
    def _meme_dir(self) -> Path:
        return self.fm.result_dir / f"{self.fm.proj}.refine.motif.meme"

    @property
    def _fimo_dir(self) -> Path:
        return self.fm.result_dir / f"{self.fm.proj}.refine.motif.fimo"
