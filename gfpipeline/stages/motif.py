"""Motif stage: discover and filter motifs using MEME/FIMO."""

from __future__ import annotations

import logging
from pathlib import Path

from gfpipeline.config.schema import PipelineConfig
from gfpipeline.core.exceptions import StageInputError
from gfpipeline.core.file_manager import FileManager
from gfpipeline.core.runner import ToolRunner
from gfpipeline.genome_db.gene_index import strip_id_prefix

log = logging.getLogger(__name__)


class MotifStage:
    """Motif discovery and genome-wide scanning stage.

    Steps:
        1. Run MEME on candidate proteins to discover motifs
        2. Run FIMO to scan genome-wide proteins with discovered motifs
        3. Parse FIMO output to build gene → motif mapping
        4. Filter genes by motif presence (any/all/min_count)
        5. Write candidates idlist and summary TSV
    """

    def __init__(self, config: PipelineConfig, runner: ToolRunner, fm: FileManager) -> None:
        self.config = config
        self.runner = runner
        self.fm = fm
        self._motif = config.motif

    # ------------------------------------------------------------------
    # Path helpers
    # ------------------------------------------------------------------

    @property
    def _pep_fa(self) -> Path:
        return self.fm.result("identify", "candidates.pep", "fa")

    @property
    def _meme_dir(self) -> Path:
        return self.fm.result_dir / f"{self.fm.proj}.motif.meme"

    @property
    def _fimo_dir(self) -> Path:
        return self.fm.result_dir / f"{self.fm.proj}.motif.fimo"

    @property
    def _candidates_idlist(self) -> Path:
        return self.fm.result("motif-filter", "candidates", "idlist")

    @property
    def _summary_tsv(self) -> Path:
        return self.fm.result("motif-filter", "summary", "tsv")

    # ------------------------------------------------------------------
    # Core methods
    # ------------------------------------------------------------------

    def run_meme(self) -> Path:
        """Call meme on candidate proteins, return path to meme.xml."""
        pep = self._pep_fa
        out_dir = self._meme_dir
        cmd = [
            self.config.tools.meme,
            str(pep),
            "-protein",
            "-nmotifs", str(self._motif.num_motifs),
            "-minw", str(self._motif.min_width),
            "-maxw", str(self._motif.max_width),
            "-oc", str(out_dir),
        ]
        log.info("Running MEME: %s", " ".join(cmd))
        self.runner.run(cmd)
        meme_xml = out_dir / "meme.xml"
        log.info("MEME complete → %s", meme_xml)
        return meme_xml

    def run_fimo(self, meme_dir: Path) -> Path:
        """Call fimo to scan genome-wide proteins, return path to fimo output dir."""
        meme_xml = meme_dir / "meme.xml"
        pep_db = self.fm.rep_pep
        out_dir = self._fimo_dir
        cmd = [
            self.config.tools.fimo,
            "--thresh", str(self._motif.fimo_pvalue),
            "--oc", str(out_dir),
            str(meme_xml),
            str(pep_db),
        ]
        log.info("Running FIMO: %s", " ".join(cmd))
        self.runner.run(cmd)
        log.info("FIMO complete → %s", out_dir)
        return out_dir

    def parse_fimo(self, fimo_dir: Path) -> dict[str, list[str]]:
        """Parse fimo.tsv, return gene_id -> [motif_id, ...] mapping."""
        fimo_tsv = fimo_dir / "fimo.tsv"
        gene_motifs: dict[str, list[str]] = {}

        text = fimo_tsv.read_text()
        for line in text.splitlines():
            line = line.strip()
            # Skip comment lines and empty lines
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 10:
                continue
            # Skip header row
            if parts[0] == "motif_id":
                continue
            motif_id = parts[0]
            sequence_name = parts[2]
            gene_id = strip_id_prefix(sequence_name)  # strip transcript:/gene: prefix if present
            if gene_id not in gene_motifs:
                gene_motifs[gene_id] = []
            if motif_id not in gene_motifs[gene_id]:
                gene_motifs[gene_id].append(motif_id)

        return gene_motifs

    def filter_genes(
        self,
        gene_motifs: dict[str, list[str]],
        all_motif_ids: list[str],
    ) -> list[str]:
        """Filter genes by motif presence.

        Modes:
            any:       gene hits any motif
            all:       gene hits all motifs in all_motif_ids
            min_count: gene hits >= min_motif_count distinct motifs
        """
        mode = self._motif.filter_mode
        result: list[str] = []

        for gene_id, motifs in gene_motifs.items():
            motif_set = set(motifs)
            if mode == "any":
                if motif_set:
                    result.append(gene_id)
            elif mode == "all":
                if set(all_motif_ids).issubset(motif_set):
                    result.append(gene_id)
            elif mode == "min_count":
                if len(motif_set) >= self._motif.min_motif_count:
                    result.append(gene_id)

        return sorted(result)

    def _parse_fimo_rows(self, fimo_dir: Path) -> list[dict]:
        """Parse fimo.tsv into list of row dicts for summary output."""
        fimo_tsv = fimo_dir / "fimo.tsv"
        rows: list[dict] = []

        text = fimo_tsv.read_text()
        for line in text.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 10:
                continue
            if parts[0] == "motif_id":
                continue
            rows.append({
                "motif_id": parts[0],
                "motif_alt_id": parts[1],
                "sequence_name": parts[2],
                "start": parts[3],
                "stop": parts[4],
                "strand": parts[5],
                "score": parts[6],
                "p_value": parts[7],
                "q_value": parts[8],
                "matched_sequence": parts[9],
            })
        return rows

    def run_filter_only(self) -> None:
        """Re-run parse + filter using existing fimo.tsv, skip MEME/FIMO."""
        fimo_dir = self._fimo_dir
        fimo_tsv = fimo_dir / "fimo.tsv"
        if not fimo_tsv.exists():
            raise StageInputError(
                f"FIMO output '{fimo_tsv}' does not exist. "
                "Please run the motif stage first."
            )
        self.fm.ensure_dirs()

        gene_motifs = self.parse_fimo(fimo_dir)
        log.info("FIMO hits: %d genes", len(gene_motifs))

        all_motif_ids: list[str] = []
        seen: set[str] = set()
        for motifs in gene_motifs.values():
            for m in motifs:
                if m not in seen:
                    all_motif_ids.append(m)
                    seen.add(m)

        filtered_ids = self.filter_genes(gene_motifs, all_motif_ids)
        log.info("Filtered %d genes by motif (%s mode)", len(filtered_ids), self._motif.filter_mode)

        self._candidates_idlist.write_text("\n".join(filtered_ids) + "\n")
        log.info("Candidates written → %s", self._candidates_idlist)

        filtered_set = set(filtered_ids)
        fimo_rows = self._parse_fimo_rows(fimo_dir)
        summary_lines = ["gene_id\tmotif_id\tmotif_name\tsequence_name\tstart\tstop\tscore\tp_value"]
        for row in fimo_rows:
            gene_id = row["sequence_name"]
            if gene_id in filtered_set:
                summary_lines.append(
                    f"{gene_id}\t{row['motif_id']}\t{row['motif_alt_id']}"
                    f"\t{row['sequence_name']}\t{row['start']}\t{row['stop']}"
                    f"\t{row['score']}\t{row['p_value']}"
                )
        self._summary_tsv.write_text("\n".join(summary_lines) + "\n")
        log.info("Summary written → %s", self._summary_tsv)
        log.info("Motif filter-only complete.")

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------

    def run(self, force: bool = False) -> None:
        """Execute the full motif stage."""
        if not self._pep_fa.exists():
            raise StageInputError(
                f"Input file '{self._pep_fa}' does not exist. "
                "Please run the identify stage first."
            )

        # Skip if outputs already exist and force=False
        if (
            self.fm.skip_if_exists(self._candidates_idlist, force)
            and self.fm.skip_if_exists(self._summary_tsv, force)
        ):
            return

        self.fm.ensure_dirs()

        # Step 1: run MEME
        meme_xml = self.run_meme()
        meme_dir = meme_xml.parent

        # Step 2: run FIMO
        fimo_dir = self.run_fimo(meme_dir)

        # Step 3: parse FIMO results
        gene_motifs = self.parse_fimo(fimo_dir)
        log.info("FIMO hits: %d genes", len(gene_motifs))

        # Collect all motif IDs discovered
        all_motif_ids: list[str] = []
        seen: set[str] = set()
        for motifs in gene_motifs.values():
            for m in motifs:
                if m not in seen:
                    all_motif_ids.append(m)
                    seen.add(m)

        # Step 4: filter genes
        filtered_ids = self.filter_genes(gene_motifs, all_motif_ids)
        log.info("Filtered %d genes by motif (%s mode)", len(filtered_ids), self._motif.filter_mode)

        # Step 5: write candidates idlist
        self._candidates_idlist.write_text("\n".join(filtered_ids) + "\n")
        log.info("Candidates written → %s", self._candidates_idlist)

        # Step 6: write summary TSV
        filtered_set = set(filtered_ids)
        fimo_rows = self._parse_fimo_rows(fimo_dir)
        summary_lines = ["gene_id\tmotif_id\tmotif_name\tsequence_name\tstart\tstop\tscore\tp_value"]
        for row in fimo_rows:
            gene_id = row["sequence_name"]
            if gene_id in filtered_set:
                summary_lines.append(
                    f"{gene_id}\t{row['motif_id']}\t{row['motif_alt_id']}"
                    f"\t{row['sequence_name']}\t{row['start']}\t{row['stop']}"
                    f"\t{row['score']}\t{row['p_value']}"
                )
        self._summary_tsv.write_text("\n".join(summary_lines) + "\n")
        log.info("Summary written → %s", self._summary_tsv)

        log.info("Motif stage complete.")
