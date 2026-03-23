"""Properties stage: calculate physicochemical properties of candidate proteins."""

from __future__ import annotations

import logging
from pathlib import Path

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from gfpipeline.config.schema import PipelineConfig
from gfpipeline.core.exceptions import StageInputError
from gfpipeline.core.file_manager import FileManager

log = logging.getLogger(__name__)

# Standard amino acid order for output columns
_AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")

_ONLINE_TOOL_HINTS = """\
理化性质计算完成。以下在线工具可用于进一步分析：
- 跨膜域预测（DeepTMHMM）：https://dtu.biolib.com/DeepTMHMM
- 信号肽预测（SignalP-6.0）：https://services.healthtech.dtu.dk/services/SignalP-6.0/
- 亚细胞定位预测（WoLF PSORT）：https://wolfpsort.hgc.jp/
输入文件：{pep_fa_path}"""

# TSV column order (before AA composition columns)
_SCALAR_COLS = [
    "gene_id",
    "length",
    "molecular_weight",
    "isoelectric_point",
    "instability_index",
    "aliphatic_index",
    "gravy",
    "aromaticity",
]


def _aliphatic_index(aa_percent: dict[str, float]) -> float:
    """Aliphatic index = X(Ala) + 2.9*X(Val) + 3.9*(X(Ile) + X(Leu)).

    X values are mole fractions (0–1), so divide percent by 100.
    Reference: Ikai (1980) J Biochem 88:1895-1898.
    """
    ala = aa_percent.get("A", 0.0) / 100.0
    val = aa_percent.get("V", 0.0) / 100.0
    ile = aa_percent.get("I", 0.0) / 100.0
    leu = aa_percent.get("L", 0.0) / 100.0
    return (ala + 2.9 * val + 3.9 * (ile + leu)) * 100.0


class PropertiesStage:
    """Calculate physicochemical properties of candidate proteins."""

    def __init__(self, config: PipelineConfig, fm: FileManager) -> None:
        self.config = config
        self.fm = fm

    # ------------------------------------------------------------------
    # Path helpers
    # ------------------------------------------------------------------

    @property
    def _pep_fa(self) -> Path:
        return self.fm.result("identify", "candidates.pep", "fa")

    @property
    def _output_tsv(self) -> Path:
        return self.fm.result("properties", "properties", "tsv")

    # ------------------------------------------------------------------
    # Core methods
    # ------------------------------------------------------------------

    def calc_properties(self, pep_fa: Path) -> list[dict]:
        """Calculate physicochemical properties for each protein sequence.

        Properties computed (all locally, no network required):
        - length: number of amino acids
        - molecular_weight: Da
        - isoelectric_point: pI
        - instability_index: Guruprasad et al. 1990
        - aliphatic_index: Ikai 1980
        - gravy: Grand Average of Hydropathicity (Kyte & Doolittle 1982)
        - aromaticity: fraction of Phe+Trp+Tyr
        - amino acid composition (%): one column per standard AA

        Args:
            pep_fa: Path to input FASTA file with protein sequences.

        Returns:
            List of dicts, one per sequence.
        """
        results: list[dict] = []
        for record in SeqIO.parse(pep_fa, "fasta"):
            seq_str = str(record.seq).replace("*", "")  # remove stop codon
            if not seq_str:
                log.warning("Skipping empty sequence: %s", record.id)
                continue
            analysis = ProteinAnalysis(seq_str)
            aa_pct = analysis.amino_acids_percent

            row: dict = {
                "gene_id": record.id,
                "length": len(seq_str),
                "molecular_weight": analysis.molecular_weight(),
                "isoelectric_point": analysis.isoelectric_point(),
                "instability_index": analysis.instability_index(),
                "aliphatic_index": _aliphatic_index(aa_pct),
                "gravy": analysis.gravy(),
                "aromaticity": analysis.aromaticity(),
            }
            for aa in _AA_ORDER:
                row[aa] = aa_pct.get(aa, 0.0)

            results.append(row)

        return results

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------

    def run(self, force: bool = False) -> None:
        """Calculate properties → properties.tsv, then print online tool hints."""
        pep_fa = self._pep_fa
        if not pep_fa.exists():
            raise StageInputError(
                f"Input file '{pep_fa}' does not exist. "
                "Please run the identify stage first."
            )

        output = self._output_tsv
        if self.fm.skip_if_exists(output, force):
            return

        self.fm.ensure_dirs()

        rows = self.calc_properties(pep_fa)
        log.info("Calculated properties for %d sequences", len(rows))

        # Write TSV
        header = "\t".join(_SCALAR_COLS + _AA_ORDER)
        lines = [header]
        for row in rows:
            values = (
                [row["gene_id"]]
                + [str(row["length"])]
                + [f"{row[k]:.4f}" for k in _SCALAR_COLS[2:]]
                + [f"{row[aa]:.6f}" for aa in _AA_ORDER]
            )
            lines.append("\t".join(values))

        output.write_text("\n".join(lines) + "\n")
        log.info("Properties written → %s", output)

        # Print online tool hints
        print(_ONLINE_TOOL_HINTS.format(pep_fa_path=pep_fa))
