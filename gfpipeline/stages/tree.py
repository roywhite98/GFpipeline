"""Tree stage: phylogenetic tree construction via muscle → trimal → iqtree2."""

from __future__ import annotations

import logging

from gfpipeline.config.schema import PipelineConfig
from gfpipeline.core.exceptions import StageInputError
from gfpipeline.core.file_manager import FileManager
from gfpipeline.core.runner import ToolRunner

log = logging.getLogger(__name__)


class TreeStage:
    """Phylogenetic tree construction stage.

    Steps:
        1. Check candidates.pep.fa exists
        2. muscle multiple sequence alignment → tree.pep.afa
        3. trimal -automated1 → tree.pep.trimed.afa
        4. iqtree2 -bb {bootstrap} -bnni -nt AUTO → tree files
    """

    def __init__(self, config: PipelineConfig, runner: ToolRunner, fm: FileManager) -> None:
        self.config = config
        self.runner = runner
        self.fm = fm
        self._tools = config.tools
        self._tree = config.tree

    # ------------------------------------------------------------------
    # Internal path helpers
    # ------------------------------------------------------------------

    @property
    def _pep_fa(self):
        return self.fm.result("identify", "candidates.pep", "fa")

    @property
    def _tree_afa(self):
        return self.fm.result("tree", "pep", "afa")

    @property
    def _tree_trimed_afa(self):
        return self.fm.result("tree", "pep.trimed", "afa")

    # ------------------------------------------------------------------
    # Public methods
    # ------------------------------------------------------------------

    def run(self, force: bool = False) -> None:
        """Execute the full tree stage."""
        # Step 1: check input
        if not self._pep_fa.exists():
            raise StageInputError(
                f"Input file '{self._pep_fa}' does not exist. "
                "Please run the identify stage first."
            )

        self.fm.ensure_dirs()

        # Step 2: muscle alignment
        if not self.fm.skip_if_exists(self._tree_afa, force):
            log.info("Running muscle alignment → %s", self._tree_afa)
            self.runner.run([
                self._tools.muscle,
                "-align", str(self._pep_fa),
                "-output", str(self._tree_afa),
            ])
            log.info("Muscle alignment done → %s", self._tree_afa)

        # Step 3: trimal
        if not self.fm.skip_if_exists(self._tree_trimed_afa, force):
            log.info("Running trimal → %s", self._tree_trimed_afa)
            self.runner.run([
                self._tools.trimal,
                "-in", str(self._tree_afa),
                "-out", str(self._tree_trimed_afa),
                "-automated1",
            ])
            log.info("trimal done → %s", self._tree_trimed_afa)

        # Step 4: iqtree2 (output prefix = trimmed afa path)
        log.info("Running iqtree2 with bootstrap=%d", self._tree.iqtree_bootstrap)
        self.runner.run([
            self._tools.iqtree,
            "-s", str(self._tree_trimed_afa),
            "-bb", str(self._tree.iqtree_bootstrap),
            "-bnni",
            "-nt", "AUTO",
        ])
        log.info("iqtree2 done. Tree files at %s.*", self._tree_trimed_afa)

        log.info("Tree stage complete.")
