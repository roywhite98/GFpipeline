"""Identify stage: gene family identification via HMM search and BLAST search."""

from __future__ import annotations

import logging
from pathlib import Path

from gfpipeline.config.schema import PipelineConfig
from gfpipeline.core.exceptions import StageInputError
from gfpipeline.core.file_manager import FileManager
from gfpipeline.core.runner import ToolRunner
from gfpipeline.core.sequence import (
    extract_fasta,
    parse_blast_idlist,
    parse_hmm_idlist,
)

log = logging.getLogger(__name__)


class IdentifyStage:
    """Gene family identification stage.

    Data flow:
        ref.fa → muscle → ref.afa → hmmbuild → {Proj}.hmm
                                                  ↓
                                            hmmsearch → hmm.out → hmm.idlist
                                                  ↓
                                            hmmemit → hmmemit.out → blastp → blast.out → blast.idlist
                                                                                  ↓
                                            merge(hmm.idlist, blast.idlist) → candidates.idlist
                                                                                  ↓
                                            transcript→gene → candidates.gene.idlist
                                                                                  ↓
                                            extract_fasta(cds/pep/genome) → candidates.{cds,pep,genome}.fa
    """

    def __init__(self, config: PipelineConfig, runner: ToolRunner, fm: FileManager) -> None:
        self.config = config
        self.runner = runner
        self.fm = fm
        self._proj = config.project_name
        self._tools = config.tools
        self._identify = config.identify
        self._db = config.databases

    # ------------------------------------------------------------------
    # Internal path helpers
    # ------------------------------------------------------------------

    @property
    def _ref_fa(self) -> Path:
        return Path(self.config.databases.ref_fa)

    @property
    def _ref_stem(self) -> str:
        """Stem of ref_fa filename, used to name derived files (hmm, hmmemit).

        e.g. wrky.ref.fa  → "wrky"
             ARF.ref.fa   → "ARF"
             myseqs.fa    → "myseqs"
        """
        name = self._ref_fa.name          # e.g. "wrky.ref.fa"
        # Strip known double-extension patterns like ".ref.fa" / ".ref.fasta"
        for suffix in (".ref.fa", ".ref.fasta", ".ref.faa"):
            if name.endswith(suffix):
                return name[: -len(suffix)]
        # Fallback: strip last extension only
        return self._ref_fa.stem

    @property
    def _hmm(self) -> Path:
        return self._ref_fa.parent / f"{self._ref_stem}.hmm"

    @property
    def _hmmemit_out(self) -> Path:
        return self._ref_fa.parent / f"{self._ref_stem}.hmmemit.out"

    @property
    def _ref_afa(self) -> Path:
        return self.fm.result("identify", "ref", "afa")

    @property
    def _hmm_out(self) -> Path:
        return self.fm.result("identify", "hmm", "out")

    @property
    def _hmm_idlist(self) -> Path:
        return self.fm.result("identify", "hmm", "idlist")

    @property
    def _blast_out(self) -> Path:
        return self.fm.result("identify", "blast", "out")

    @property
    def _blast_idlist(self) -> Path:
        return self.fm.result("identify", "blast", "idlist")

    @property
    def _candidates_idlist(self) -> Path:
        return self.fm.result("identify", "candidates", "idlist")

    @property
    def _candidates_gene_idlist(self) -> Path:
        return self.fm.result("identify", "candidates.gene", "idlist")

    @property
    def _candidates_cds_fa(self) -> Path:
        return self.fm.result("identify", "candidates.cds", "fa")

    @property
    def _candidates_pep_fa(self) -> Path:
        return self.fm.result("identify", "candidates.pep", "fa")

    @property
    def _candidates_genome_fa(self) -> Path:
        return self.fm.result("identify", "candidates.genome", "fa")

    # ------------------------------------------------------------------
    # Public methods
    # ------------------------------------------------------------------

    def build_hmm(self) -> None:
        """muscle multiple sequence alignment → hmmbuild HMM profile."""
        # Step 1: muscle alignment
        if not self.fm.skip_if_exists(self._ref_afa, force=False):
            log.info("Running muscle alignment on %s", self._ref_fa)
            self.runner.run([
                self._tools.muscle,
                "-align", str(self._ref_fa),
                "-output", str(self._ref_afa),
            ])
            log.info("Muscle alignment done → %s", self._ref_afa)

        # Step 2: hmmbuild
        if not self.fm.skip_if_exists(self._hmm, force=False):
            log.info("Running hmmbuild → %s", self._hmm)
            self.runner.run([
                self._tools.hmmbuild,
                str(self._hmm),
                str(self._ref_afa),
            ])
            log.info("hmmbuild done → %s", self._hmm)

    def hmm_search(self) -> list[str]:
        """hmmsearch → parse results → return candidate ID list."""
        if not self.fm.skip_if_exists(self._hmm_out, force=False):
            log.info("Running hmmsearch with %s", self._hmm)
            self.runner.run([
                self._tools.hmmsearch,
                "--noali",
                "-o", str(self._hmm_out),
                str(self._hmm),
                str(self.fm.rep_pep),
            ])
            log.info("hmmsearch done → %s", self._hmm_out)

        if not self.fm.skip_if_exists(self._hmm_idlist, force=False):
            ids = parse_hmm_idlist(self._hmm_out, self._identify.hmm_evalue)
            self._hmm_idlist.write_text("\n".join(ids) + ("\n" if ids else ""))
            log.info("HMM candidates: %d IDs → %s", len(ids), self._hmm_idlist)
        else:
            ids = [
                line.strip()
                for line in self._hmm_idlist.read_text().splitlines()
                if line.strip()
            ]

        return ids

    def blast_search(self) -> list[str]:
        """hmmemit → blastp → parse tabular output → return candidate ID list."""
        # Step 1: hmmemit
        if not self.fm.skip_if_exists(self._hmmemit_out, force=False):
            log.info("Running hmmemit → %s", self._hmmemit_out)
            self.runner.run([
                self._tools.hmmemit,
                "-o", str(self._hmmemit_out),
                str(self._hmm),
            ])
            log.info("hmmemit done → %s", self._hmmemit_out)

        # Step 2: blastp
        if not self.fm.skip_if_exists(self._blast_out, force=False):
            log.info("Running blastp → %s", self._blast_out)
            self.runner.run([
                self._tools.blastp,
                "-query", str(self._hmmemit_out),
                "-db", self.fm.blast_db_prefix,
                "-outfmt", "6",
                "-num_threads", str(self._identify.blast_threads),
                "-evalue", str(self._identify.blast_evalue),
                "-out", str(self._blast_out),
            ])
            log.info("blastp done → %s", self._blast_out)

        # Step 3: parse
        if not self.fm.skip_if_exists(self._blast_idlist, force=False):
            ids = parse_blast_idlist(self._blast_out)
            self._blast_idlist.write_text("\n".join(ids) + ("\n" if ids else ""))
            log.info("BLAST candidates: %d IDs → %s", len(ids), self._blast_idlist)
        else:
            ids = [
                line.strip()
                for line in self._blast_idlist.read_text().splitlines()
                if line.strip()
            ]

        return ids

    def merge_ids(self, hmm_ids: list[str], blast_ids: list[str]) -> list[str]:
        """Merge, sort, and deduplicate two ID lists."""
        merged = sorted(set(hmm_ids) | set(blast_ids))
        return merged

    def extract_sequences(self, ids: list[str], gene_ids: list[str]) -> None:
        """Extract protein sequences for the candidate IDs (tree/domain/motif stages need this)."""
        # Protein sequences (by transcript ID) - required by tree/domain/motif stages
        if not self.fm.skip_if_exists(self._candidates_pep_fa, force=False):
            count = extract_fasta(self.fm.rep_pep, ids, self._candidates_pep_fa)
            log.info("Extracted %d protein sequences → %s", count, self._candidates_pep_fa)

    def run(self, force: bool = False) -> None:
        """Execute the full identify stage."""
        # Validate inputs
        if not self._ref_fa.exists() and not self._hmm.exists():
            raise StageInputError(
                f"Neither reference FASTA '{self._ref_fa}' nor HMM profile '{self._hmm}' "
                "exists. Please provide at least one of these inputs."
            )

        if not self.fm.rep_pep.exists():
            raise StageInputError(
                f"Representative protein FASTA not found: {self.fm.rep_pep}. "
                "Please run 'gfpipeline genome-db' first."
            )

        self.fm.ensure_dirs()

        # Build HMM if ref.fa exists and hmm doesn't
        if self._ref_fa.exists() and not self._hmm.exists():
            self.build_hmm()

        # HMM search
        hmm_ids = self.hmm_search()

        # BLAST search
        blast_ids = self.blast_search()

        # Merge IDs
        merged_ids = self.merge_ids(hmm_ids, blast_ids)
        if not self.fm.skip_if_exists(self._candidates_idlist, force):
            self._candidates_idlist.write_text(
                "\n".join(merged_ids) + ("\n" if merged_ids else "")
            )
            log.info("Merged candidates: %d IDs → %s", len(merged_ids), self._candidates_idlist)

        # Convert transcript IDs to gene IDs via gene2transcript.tsv reverse mapping
        gene2transcript_tsv = Path(self.config.genome_db.index_dir) / "gene2transcript.tsv"
        transcript_to_gene: dict[str, str] = {}
        if gene2transcript_tsv.exists():
            for line in gene2transcript_tsv.read_text().splitlines():
                if line.startswith("gene_id") or not line.strip():
                    continue
                parts_line = line.split("\t")
                if len(parts_line) >= 2:
                    transcript_to_gene[parts_line[1]] = parts_line[0]

        gene_ids = sorted(set(
            transcript_to_gene.get(t, t) for t in merged_ids
        ))
        if not self.fm.skip_if_exists(self._candidates_gene_idlist, force):
            self._candidates_gene_idlist.write_text(
                "\n".join(gene_ids) + ("\n" if gene_ids else "")
            )
            log.info(
                "Gene IDs: %d → %s", len(gene_ids), self._candidates_gene_idlist
            )

        # Extract sequences
        self.extract_sequences(merged_ids, gene_ids)

        log.info("Identify stage complete.")
