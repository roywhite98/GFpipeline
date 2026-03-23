"""Collinearity stage: all-vs-all BLAST, synteny analysis, block extraction."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from gfpipeline.config.schema import PipelineConfig
from gfpipeline.core.exceptions import StageInputError
from gfpipeline.core.file_manager import FileManager
from gfpipeline.core.runner import ToolRunner

log = logging.getLogger(__name__)


class CollinearityStage:
    """Collinearity / synteny analysis stage.

    Steps:
        1. All-vs-all blastp self-comparison
        2. Run JCVI or MCScanX for synteny detection
        3. Extract collinearity blocks containing target genes
        4. Parse GFF3 for gene locations
    """

    def __init__(self, config: PipelineConfig, runner: ToolRunner, fm: FileManager) -> None:
        self.config = config
        self.runner = runner
        self.fm = fm

    # ------------------------------------------------------------------
    # Path helpers
    # ------------------------------------------------------------------

    @property
    def _gene_idlist(self) -> Path:
        return self.fm.result("identify", "candidates.gene", "idlist")

    @property
    def _blast_out(self) -> Path:
        return self.fm.result("collinearity", "blast", "out")

    @property
    def _collinearity_dir(self) -> Path:
        return self.fm.result_dir / f"{self.fm.proj}.collinearity"

    @property
    def _blocks_tsv(self) -> Path:
        return self.fm.result("collinearity", "blocks", "tsv")

    @property
    def _gene_location_tsv(self) -> Path:
        return self.fm.result("collinearity", "gene-location", "tsv")

    @property
    def _genome_name(self) -> str:
        return self.config.genome_db.genome_name or self.config.project_name

    # ------------------------------------------------------------------
    # Input validation
    # ------------------------------------------------------------------

    def _validate_inputs(self) -> None:
        """Raise StageInputError if any required input file is missing."""
        missing = []
        for path, label in [
            (self._gene_idlist, "candidates.gene.idlist"),
            (Path(self.config.databases.gff3), "databases.gff3"),
            (self.fm.rep_pep, "rep_pep (run genome-db first)"),
        ]:
            if not path.exists():
                missing.append(f"'{path}' ({label})")
        if missing:
            raise StageInputError(
                "Required input file(s) missing: " + ", ".join(missing)
            )

    # ------------------------------------------------------------------
    # Core methods
    # ------------------------------------------------------------------

    def run_all_vs_all_blast(self) -> Path:
        """All-vs-all blastp self-comparison, output collinearity.blast.out."""
        pep = self.fm.rep_pep
        blast_db = self.fm.blast_db_prefix
        threads = self.config.collinearity.blast_threads
        blast_out = self._blast_out

        cmd = [
            self.config.tools.blastp,
            "-query", str(pep),
            "-db", str(blast_db),
            "-outfmt", "6",
            "-num_threads", str(threads),
            "-out", str(blast_out),
        ]
        log.info("Running all-vs-all blastp: %s", " ".join(cmd))
        self.runner.run(cmd)
        log.info("BLAST complete → %s", blast_out)
        return blast_out

    def run_jcvi(self) -> None:
        """Call python -m jcvi.compara.catalog ortholog, output to collinearity/ dir."""
        genome_name = self._genome_name
        threads = self.config.collinearity.blast_threads
        out_dir = self._collinearity_dir
        out_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            sys.executable,
            "-m", "jcvi.compara.catalog",
            "ortholog",
            genome_name,
            genome_name,
            f"--cpu={threads}",
        ]
        log.info("Running JCVI: %s", " ".join(cmd))
        self.runner.run(cmd, cwd=str(out_dir))
        log.info("JCVI complete → %s", out_dir)

    def run_mcscanx(self) -> None:
        """Call MCScanX, output to collinearity/ dir."""
        genome_name = self._genome_name
        out_dir = self._collinearity_dir
        out_dir.mkdir(parents=True, exist_ok=True)

        input_prefix = str(out_dir / genome_name)
        cmd = [
            self.config.tools.mcscanx,
            input_prefix,
        ]
        log.info("Running MCScanX: %s", " ".join(cmd))
        self.runner.run(cmd, cwd=str(out_dir))
        log.info("MCScanX complete → %s", out_dir)

    def extract_target_blocks(self, gene_ids: list[str]) -> list[dict]:
        """Extract collinearity blocks containing target genes.

        Reads any .collinearity or .synteny files in the collinearity dir.
        Writes collinearity.blocks.tsv with columns:
            block_id, gene1, gene2, chromosome1, chromosome2, score
        Returns list of block dicts.
        """
        gene_id_set = set(gene_ids)
        blocks: list[dict] = []

        col_dir = self._collinearity_dir
        col_files: list[Path] = []
        if col_dir.exists():
            col_files = list(col_dir.glob("*.collinearity")) + list(col_dir.glob("*.synteny"))

        if col_files:
            current_block_id: str | None = None
            current_score: str = "."
            for col_file in col_files:
                for line in col_file.read_text().splitlines():
                    line = line.strip()
                    if not line or line.startswith("#"):
                        # Try to parse block header for score
                        if line.startswith("## Alignment"):
                            parts = line.split()
                            # e.g. ## Alignment 0: score=...
                            for part in parts:
                                if part.startswith("score="):
                                    current_score = part.split("=", 1)[1].rstrip(",")
                            # Extract block id
                            for part in parts:
                                if part.endswith(":"):
                                    current_block_id = part.rstrip(":")
                        continue
                    parts = line.split()
                    if len(parts) < 3:
                        continue
                    # Lines like: index  gene1  gene2  [e-value]
                    gene1 = parts[1] if len(parts) > 1 else ""
                    gene2 = parts[2] if len(parts) > 2 else ""
                    if gene1 in gene_id_set or gene2 in gene_id_set:
                        # Derive chromosome from gene id (prefix before last dot or underscore)
                        chr1 = gene1.rsplit(".", 1)[0] if "." in gene1 else gene1
                        chr2 = gene2.rsplit(".", 1)[0] if "." in gene2 else gene2
                        blocks.append({
                            "block_id": current_block_id or ".",
                            "gene1": gene1,
                            "gene2": gene2,
                            "chromosome1": chr1,
                            "chromosome2": chr2,
                            "score": current_score,
                        })

        # Write TSV
        header = "block_id\tgene1\tgene2\tchromosome1\tchromosome2\tscore"
        lines = [header]
        for b in blocks:
            lines.append(
                f"{b['block_id']}\t{b['gene1']}\t{b['gene2']}"
                f"\t{b['chromosome1']}\t{b['chromosome2']}\t{b['score']}"
            )
        self._blocks_tsv.write_text("\n".join(lines) + "\n")
        log.info("Blocks written → %s (%d rows)", self._blocks_tsv, len(blocks))
        return blocks

    def parse_gene_locations(self, gff3: Path, gene_ids: list[str]) -> None:
        """Parse GFF3 for gene features matching gene_ids.

        Writes collinearity.gene-location.tsv with columns:
            gene_id, chromosome, start, end, strand
        """
        gene_id_set = set(gene_ids)
        rows: list[dict] = []

        for line in gff3.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            if feature_type != "gene":
                continue
            chrom = parts[0]
            start = parts[3]
            end = parts[4]
            strand = parts[6]
            attrs = parts[8]

            # Parse gene ID from attributes (ID=xxx or Name=xxx)
            gene_id: str | None = None
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("ID="):
                    gene_id = attr[3:]
                    break
            if gene_id is None:
                continue

            if gene_id in gene_id_set:
                rows.append({
                    "gene_id": gene_id,
                    "chromosome": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                })

        header = "gene_id\tchromosome\tstart\tend\tstrand"
        lines = [header]
        for r in rows:
            lines.append(
                f"{r['gene_id']}\t{r['chromosome']}\t{r['start']}\t{r['end']}\t{r['strand']}"
            )
        self._gene_location_tsv.write_text("\n".join(lines) + "\n")
        log.info("Gene locations written → %s (%d genes)", self._gene_location_tsv, len(rows))

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------

    def run(self, force: bool = False) -> None:
        """Execute the full collinearity stage."""
        self._validate_inputs()

        # Skip if outputs already exist and force=False
        if (
            self.fm.skip_if_exists(self._blocks_tsv, force)
            and self.fm.skip_if_exists(self._gene_location_tsv, force)
        ):
            return

        self.fm.ensure_dirs()
        self._collinearity_dir.mkdir(parents=True, exist_ok=True)

        # Step 1: all-vs-all BLAST
        self.run_all_vs_all_blast()

        # Step 2: synteny tool
        tool = self.config.collinearity.tool
        if tool == "jcvi":
            self.run_jcvi()
        else:
            self.run_mcscanx()

        # Step 3: read gene ids
        gene_ids = [
            g.strip()
            for g in self._gene_idlist.read_text().splitlines()
            if g.strip()
        ]

        # Step 4: extract blocks
        self.extract_target_blocks(gene_ids)

        # Step 5: parse gene locations
        self.parse_gene_locations(Path(self.config.databases.gff3), gene_ids)

        log.info("Collinearity stage complete.")
