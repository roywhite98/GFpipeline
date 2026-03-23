"""Trans stage: expression matrix analysis, heatmap, DEG filtering, and Venn diagram."""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd

from gfpipeline.config.schema import PipelineConfig
from gfpipeline.core.exceptions import StageInputError
from gfpipeline.core.file_manager import FileManager

log = logging.getLogger(__name__)

try:
    from matplotlib_venn import venn2, venn3
    _VENN_AVAILABLE = True
except ImportError:
    log.warning("matplotlib_venn is not installed; plot_venn will be skipped.")
    _VENN_AVAILABLE = False


class TransStage:
    """Expression analysis stage: heatmap, DEG filtering, and Venn diagram."""

    def __init__(self, config: PipelineConfig, fm: FileManager) -> None:
        self.config = config
        self.fm = fm

    # ------------------------------------------------------------------
    # Path helpers
    # ------------------------------------------------------------------

    @property
    def _gene_idlist(self) -> Path:
        return self.fm.result("identify", "candidates.gene", "idlist")

    @property
    def _heatmap_pdf(self) -> Path:
        return self.fm.result("trans", "expression-heatmap", "pdf")

    @property
    def _deg_tsv(self) -> Path:
        return self.fm.result("trans", "deg", "tsv")

    @property
    def _venn_pdf(self) -> Path:
        return self.fm.result("trans", "venn", "pdf")

    # ------------------------------------------------------------------
    # Core methods
    # ------------------------------------------------------------------

    def load_expression_matrix(self, path: Path) -> pd.DataFrame:
        """Load expression matrix (rows=genes, cols=samples) from TSV/CSV."""
        sep = "\t" if path.suffix.lower() in (".tsv", ".txt") else ","
        df = pd.read_csv(path, sep=sep, index_col=0)
        return df

    def filter_family_members(self, matrix: pd.DataFrame, gene_ids: list[str]) -> pd.DataFrame:
        """Filter rows to only include gene family members."""
        mask = matrix.index.isin(gene_ids)
        return matrix.loc[mask]

    def plot_heatmap(self, matrix: pd.DataFrame, output: Path) -> None:
        """Draw seaborn clustermap, save to PDF."""
        g = sns.clustermap(matrix, cmap="viridis", figsize=(10, 8))
        g.savefig(output)
        plt.close("all")
        log.info("Heatmap written → %s", output)

    def filter_deg(self, matrix: pd.DataFrame) -> pd.DataFrame:
        """Filter by logfc_threshold and pvalue_threshold from config."""
        logfc_threshold = self.config.trans.logfc_threshold
        pvalue_threshold = self.config.trans.pvalue_threshold

        if logfc_threshold is None and pvalue_threshold is None:
            return matrix

        result = matrix.copy()

        if logfc_threshold is not None:
            # Look for a dedicated logfc column first
            logfc_cols = [c for c in result.columns if c.lower() in ("logfc", "log2fc")]
            if logfc_cols:
                col = logfc_cols[0]
                result = result[result[col].abs() >= logfc_threshold]
            else:
                # Filter rows where any numeric column value >= threshold
                numeric = result.select_dtypes(include="number")
                if numeric.empty:
                    log.warning(
                        "No numeric columns found for logfc filtering; returning full matrix."
                    )
                else:
                    mask = (numeric.abs() >= logfc_threshold).any(axis=1)
                    result = result.loc[mask]

        if pvalue_threshold is not None:
            pval_cols = [c for c in result.columns if c.lower() in ("pvalue", "p_value", "padj")]
            if pval_cols:
                col = pval_cols[0]
                result = result[result[col] < pvalue_threshold]
            else:
                log.warning(
                    "No pvalue/p_value/padj column found for pvalue filtering; "
                    "returning current matrix without pvalue filter."
                )

        return result

    def plot_venn(self, groups: list[str], matrix: pd.DataFrame, output: Path) -> None:
        """Draw Venn diagram using matplotlib_venn, save to PDF."""
        if not _VENN_AVAILABLE:
            log.warning("matplotlib_venn not available; skipping Venn diagram.")
            return

        if len(groups) not in (2, 3):
            log.warning(
                "plot_venn requires exactly 2 or 3 groups, got %d; skipping.", len(groups)
            )
            return

        # Build sets: genes with expression > 0 in each group
        sets = []
        for grp in groups:
            if grp not in matrix.columns:
                log.warning("Group column '%s' not found in matrix; skipping Venn.", grp)
                return
            expressed = set(matrix.index[matrix[grp] > 0])
            sets.append(expressed)

        fig, ax = plt.subplots()
        if len(groups) == 2:
            venn2(sets, set_labels=groups, ax=ax)
        else:
            venn3(sets, set_labels=groups, ax=ax)

        fig.savefig(output)
        plt.close("all")
        log.info("Venn diagram written → %s", output)

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------

    def run(self, force: bool = False) -> None:
        """Orchestrate the full trans stage."""
        expr_path = self.config.trans.expression_matrix
        if expr_path is None:
            raise StageInputError(
                "config.trans.expression_matrix is not set. "
                "Please provide an expression matrix file."
            )

        expr_file = Path(expr_path)
        if not expr_file.exists():
            raise StageInputError(
                f"Expression matrix file '{expr_file}' does not exist."
            )

        # Check if all outputs already exist
        outputs = [self._heatmap_pdf, self._deg_tsv, self._venn_pdf]
        if all(p.exists() for p in outputs) and not force:
            log.info("All trans outputs exist; skipping (use force=True to rerun).")
            return

        self.fm.ensure_dirs()

        # Load expression matrix
        matrix = self.load_expression_matrix(expr_file)
        log.info("Loaded expression matrix: %d genes × %d samples", *matrix.shape)

        # Load gene family IDs
        gene_ids: list[str] = []
        if self._gene_idlist.exists():
            gene_ids = [
                line.strip()
                for line in self._gene_idlist.read_text().splitlines()
                if line.strip()
            ]

        # Filter to family members
        family_matrix = self.filter_family_members(matrix, gene_ids)
        log.info("Family members in matrix: %d", len(family_matrix))

        # Heatmap
        if family_matrix.empty:
            log.warning(
                "No gene family members found in expression matrix; skipping heatmap."
            )
        else:
            self.plot_heatmap(family_matrix, self._heatmap_pdf)

        # DEG filtering
        deg = self.filter_deg(matrix)
        deg.to_csv(self._deg_tsv, sep="\t")
        log.info("DEG table written → %s (%d rows)", self._deg_tsv, len(deg))

        # Venn diagram
        groups = self.config.trans.venn_groups
        if groups:
            self.plot_venn(groups, matrix, self._venn_pdf)
        else:
            log.info("No venn_groups configured; skipping Venn diagram.")
