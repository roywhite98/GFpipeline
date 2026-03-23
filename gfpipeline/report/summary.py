"""Summary reporter for gene-family-pipeline."""

from __future__ import annotations

import glob
from datetime import datetime
from pathlib import Path

from gfpipeline.config.schema import PipelineConfig
from gfpipeline.core.file_manager import FileManager

_NOT_EXECUTED = "未执行"


def _count_idlist(path: Path) -> int:
    """Count non-empty lines in an idlist file."""
    return sum(1 for line in path.read_text(encoding="utf-8").splitlines() if line.strip())


def _count_tsv_rows(path: Path) -> int:
    """Count data rows in a TSV file (total lines minus 1 header)."""
    lines = [l for l in path.read_text(encoding="utf-8").splitlines() if l.strip()]
    return max(0, len(lines) - 1)


class Summary_Reporter:
    """Collect and write a human-readable summary of pipeline outputs."""

    def __init__(self, config: PipelineConfig, fm: FileManager):
        self.config = config
        self.fm = fm

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _file_info(self, path: Path, count: int | None = None) -> str:
        """Return formatted info string for a file path."""
        if path.exists():
            if count is not None:
                return f"{path} ({count} 个基因)"
            return str(path)
        return _NOT_EXECUTED

    def _tsv_info(self, path: Path, label: str = "个块") -> str:
        """Return formatted info string for a TSV file with row count."""
        if path.exists():
            n = _count_tsv_rows(path)
            return f"{path} ({n} {label})"
        return _NOT_EXECUTED

    def _idlist_info(self, path: Path) -> str:
        """Return formatted info string for an idlist file with gene count."""
        if path.exists():
            n = _count_idlist(path)
            return f"{path} ({n} 个基因)"
        return _NOT_EXECUTED

    def _tree_info(self) -> str:
        """Find tree file (treefile or contree) and return info."""
        base = self.fm.result("tree", "pep.trimed.afa", "treefile")
        # Check exact path first
        if base.exists():
            return str(base)
        # Glob for any matching treefile or contree
        pattern_treefile = str(self.fm.result_dir / f"{self.fm.proj}.tree.pep.trimed.afa.*treefile")
        pattern_contree = str(self.fm.result_dir / f"{self.fm.proj}.tree.pep.trimed.afa.*.contree")
        for pattern in (pattern_treefile, pattern_contree):
            matches = glob.glob(pattern)
            if matches:
                return str(matches[0])
        return _NOT_EXECUTED

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def collect(self) -> dict:
        """Collect stage outputs and statistics. Returns a dict of stage -> info."""
        fm = self.fm

        # identify
        id_idlist = fm.result("identify", "candidates.gene", "idlist")
        id_pep = fm.result("identify", "candidates.pep", "fa")
        id_cds = fm.result("identify", "candidates.cds", "fa")

        # tree
        tree_file = self._tree_info()

        # domain
        domain_cdd = fm.result("domain", "cdd", "txt")

        # domain-filter
        df_idlist = fm.result("domain-filter", "candidates", "idlist")

        # motif
        motif_dir = fm.result_dir / f"{fm.proj}.motif.meme"
        mf_idlist = fm.result("motif-filter", "candidates", "idlist")

        # collinearity
        col_blocks = fm.result("collinearity", "blocks", "tsv")
        col_loc = fm.result("collinearity", "gene-location", "tsv")

        # properties
        prop_tsv = fm.result("properties", "properties", "tsv")

        # trans
        trans_heatmap = fm.result("trans", "expression-heatmap", "pdf")
        trans_deg = fm.result("trans", "deg", "tsv")

        return {
            "identify": {
                "candidates_idlist": id_idlist,
                "candidates_pep": id_pep,
                "candidates_cds": id_cds,
            },
            "tree": {
                "tree_file": tree_file,
            },
            "domain": {
                "cdd": domain_cdd,
            },
            "domain-filter": {
                "candidates_idlist": df_idlist,
            },
            "motif": {
                "meme_dir": motif_dir,
                "candidates_idlist": mf_idlist,
            },
            "collinearity": {
                "blocks_tsv": col_blocks,
                "gene_location_tsv": col_loc,
            },
            "properties": {
                "properties_tsv": prop_tsv,
            },
            "trans": {
                "heatmap_pdf": trans_heatmap,
                "deg_tsv": trans_deg,
            },
        }

    def write(self) -> Path:
        """Write summary.txt and return its path."""
        fm = self.fm
        out_path = fm.result_dir / f"{fm.proj}.summary.txt"
        data = self.collect()

        lines: list[str] = []
        lines.append("=== gene-family-pipeline 分析汇总报告 ===")
        lines.append(f"项目名称: {self.config.project_name}")
        lines.append(f"结果目录: {fm.result_dir}")
        lines.append(f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        lines.append("")

        # --- identify ---
        lines.append("--- identify 阶段 ---")
        id_idlist: Path = data["identify"]["candidates_idlist"]
        lines.append(f"候选基因 ID 列表: {self._idlist_info(id_idlist)}")
        id_pep: Path = data["identify"]["candidates_pep"]
        lines.append(f"候选蛋白序列: {self._file_info(id_pep)}")
        id_cds: Path = data["identify"]["candidates_cds"]
        lines.append(f"候选 CDS 序列: {self._file_info(id_cds)}")
        lines.append("")

        # --- tree ---
        lines.append("--- tree 阶段 ---")
        tree_info = data["tree"]["tree_file"]
        lines.append(f"系统发育树文件: {tree_info}")
        lines.append("")

        # --- domain ---
        lines.append("--- domain 阶段 ---")
        domain_cdd: Path = data["domain"]["cdd"]
        lines.append(f"CDD 结果文件: {self._file_info(domain_cdd)}")
        lines.append("")

        # --- domain-filter ---
        lines.append("--- domain-filter 阶段 ---")
        df_idlist: Path = data["domain-filter"]["candidates_idlist"]
        lines.append(f"筛选 ID 列表: {self._idlist_info(df_idlist)}")
        lines.append("")

        # --- motif ---
        lines.append("--- motif 阶段 ---")
        motif_dir: Path = data["motif"]["meme_dir"]
        lines.append(f"MEME 输出目录: {motif_dir if motif_dir.exists() else _NOT_EXECUTED}")
        mf_idlist: Path = data["motif"]["candidates_idlist"]
        lines.append(f"筛选 ID 列表: {self._idlist_info(mf_idlist)}")
        lines.append("")

        # --- collinearity ---
        lines.append("--- collinearity 阶段 ---")
        col_blocks: Path = data["collinearity"]["blocks_tsv"]
        lines.append(f"共线性块文件: {self._tsv_info(col_blocks, '个块')}")
        col_loc: Path = data["collinearity"]["gene_location_tsv"]
        lines.append(f"基因位置文件: {self._file_info(col_loc)}")
        lines.append("")

        # --- properties ---
        lines.append("--- properties 阶段 ---")
        prop_tsv: Path = data["properties"]["properties_tsv"]
        lines.append(f"理化性质表: {self._file_info(prop_tsv)}")
        lines.append("")

        # --- trans ---
        lines.append("--- trans 阶段 ---")
        trans_heatmap: Path = data["trans"]["heatmap_pdf"]
        lines.append(f"表达量热图: {self._file_info(trans_heatmap)}")
        trans_deg: Path = data["trans"]["deg_tsv"]
        lines.append(f"差异基因列表: {self._file_info(trans_deg)}")
        lines.append("")

        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text("\n".join(lines), encoding="utf-8")
        return out_path
