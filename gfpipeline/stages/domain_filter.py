"""Domain-filter stage: filter genome proteins by conserved domain."""

from __future__ import annotations

import io
import logging
import time
from dataclasses import dataclass
from pathlib import Path

import requests
from Bio import SeqIO

from gfpipeline.config.schema import PipelineConfig
from gfpipeline.core.exceptions import ApiError, StageInputError
from gfpipeline.core.file_manager import FileManager

log = logging.getLogger(__name__)

_BWRPSB_URL = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"
_BATCH_SIZE = 4000
_POLL_INTERVAL = 5
_BATCH_DELAY = 2


@dataclass
class DomainSummaryRow:
    gene_id: str
    domain_accession: str
    domain_name: str
    superfamily_accession: str
    evalue: float


class DomainFilterStage:
    """Filter genome proteins by conserved domain/superfamily.

    Steps:
        1. Check domain.cdd.txt exists
        2. Determine target domains (from config or auto-extract)
        3. Get genome CDD result (from config or submit batch)
        4. Filter genes by domain/superfamily match
        5. Output candidates.idlist and summary.tsv
    """

    def __init__(self, config: PipelineConfig, fm: FileManager) -> None:
        self.config = config
        self.fm = fm
        self._domain = config.domain

    # ------------------------------------------------------------------
    # Path helpers
    # ------------------------------------------------------------------

    @property
    def _cdd_result(self) -> Path:
        return self.fm.result("domain", "cdd", "txt")

    @property
    def _genome_cdd(self) -> Path:
        return self.fm.result("domain", "genome.cdd", "txt")

    @property
    def _candidates_idlist(self) -> Path:
        return self.fm.result("domain-filter", "candidates", "idlist")

    @property
    def _summary_tsv(self) -> Path:
        return self.fm.result("domain-filter", "summary", "tsv")

    # ------------------------------------------------------------------
    # Core methods
    # ------------------------------------------------------------------

    def extract_target_domains(self, cdd_result: Path) -> tuple[set[str], set[str]]:
        """Extract Domain_Accession and Superfamily_Accession sets from CDD result."""
        rows = self.parse_cdd_result(cdd_result)
        domains: set[str] = set()
        superfamilies: set[str] = set()
        for row in rows:
            if row.domain_accession:
                domains.add(row.domain_accession)
            if row.superfamily_accession and row.superfamily_accession != "N/A":
                superfamilies.add(row.superfamily_accession)
        return domains, superfamilies

    def parse_cdd_result(self, cdd_path: Path) -> list[DomainSummaryRow]:
        """Parse CDD result TSV file, return list of DomainSummaryRow.

        CDD TSV columns (0-indexed):
          0: Query, 1: Hit_type, 2: PSSM_ID, 3: From, 4: To,
          5: E-Value, 6: Bitscore, 7: Accession, 8: Short_name,
          9: Incomplete, 10: Superfamily
        """
        rows: list[DomainSummaryRow] = []
        try:
            text = cdd_path.read_text()
        except OSError as e:
            raise StageInputError(f"Failed to read CDD result '{cdd_path}': {e}") from e

        for line in text.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 11:
                continue
            # Skip the column header row
            if parts[0].strip() == "Query":
                continue
            query = parts[0].strip()
            # gene_id is the query ID (may have extra description after space)
            gene_id = query.split()[0] if query else ""
            evalue_str = parts[5].strip()
            accession = parts[7].strip()
            short_name = parts[8].strip()
            superfamily = parts[10].strip()
            try:
                evalue = float(evalue_str)
            except ValueError:
                evalue = 0.0
            rows.append(DomainSummaryRow(
                gene_id=gene_id,
                domain_accession=accession,
                domain_name=short_name,
                superfamily_accession=superfamily,
                evalue=evalue,
            ))
        return rows

    def filter_by_domains(
        self,
        rows: list[DomainSummaryRow],
        target_domains: set[str],
        target_superfamilies: set[str],
    ) -> list[str]:
        """Return sorted list of gene IDs that match target domains or superfamilies."""
        result: set[str] = set()
        for row in rows:
            if row.domain_accession in target_domains:
                result.add(row.gene_id)
            elif row.superfamily_accession in target_superfamilies:
                result.add(row.gene_id)
        return sorted(result)

    def submit_genome_batch(self, pep_db: Path) -> Path:
        """Submit full genome proteins in batches (≤4000), merge results to genome.cdd.txt."""
        records = list(SeqIO.parse(str(pep_db), "fasta"))
        log.info("Genome proteins: %d sequences", len(records))

        all_lines: list[str] = []
        header_written = False

        for i in range(0, len(records), _BATCH_SIZE):
            batch = records[i: i + _BATCH_SIZE]
            batch_num = i // _BATCH_SIZE + 1
            log.info("Submitting genome batch %d (%d sequences)", batch_num, len(batch))

            buf = io.StringIO()
            SeqIO.write(batch, buf, "fasta")
            batch_fasta = buf.getvalue()

            data = {
                "useid1": "true",
                "maxhit": str(self._domain.maxhit),
                "filter": "true",
                "db": "cdd",
                "evalue": str(self._domain.evalue),
                "tdata": "hits",
                "queries": batch_fasta,
            }
            resp = requests.post(_BWRPSB_URL, data=data)
            if resp.status_code != 200:
                raise ApiError(resp.status_code, resp.reason)

            rid = self._parse_rid(resp.text)
            self._poll_status(rid)

            resp = requests.post(_BWRPSB_URL, data={"tdata": "hits", "cdsid": rid})
            if resp.status_code != 200:
                raise ApiError(resp.status_code, resp.reason)

            for line in resp.text.splitlines():
                if line.startswith("#"):
                    if not header_written:
                        all_lines.append(line)
                else:
                    all_lines.append(line)
            header_written = True

            if i + _BATCH_SIZE < len(records):
                time.sleep(_BATCH_DELAY)

        self._genome_cdd.write_text("\n".join(all_lines) + "\n")
        log.info("Genome CDD result saved → %s", self._genome_cdd)
        return self._genome_cdd

    # ------------------------------------------------------------------
    # API helpers (shared with DomainStage logic)
    # ------------------------------------------------------------------

    def _parse_rid(self, text: str) -> str:
        for line in text.splitlines():
            line = line.strip()
            if line.startswith("#cdsid"):
                parts = line.split()
                if len(parts) >= 2:
                    return parts[1]
        raise ApiError(reason=f"Could not parse RID from response: {text[:200]}")

    def _poll_status(self, rid: str) -> None:
        while True:
            print("等待结果中...")
            resp = requests.post(_BWRPSB_URL, data={"tdata": "hits", "cdsid": rid})
            if resp.status_code != 200:
                raise ApiError(resp.status_code, resp.reason)
            status = self._parse_status(resp.text)
            if status == 0:
                return
            elif status == 3:
                time.sleep(_POLL_INTERVAL)
            elif status in (1, 2, 4, 5):
                raise ApiError(reason=f"CD-Search job {rid} failed with status code {status}")
            else:
                time.sleep(_POLL_INTERVAL)

    def _parse_status(self, text: str) -> int:
        for line in text.splitlines():
            line = line.strip()
            if line.startswith("#status"):
                parts = line.split()
                if len(parts) >= 2:
                    return int(parts[1])
        return 3

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------

    def run(self, force: bool = False) -> None:
        """Execute the full domain-filter stage."""
        if not self._cdd_result.exists():
            raise StageInputError(
                f"Input file '{self._cdd_result}' does not exist. "
                "Please run the domain stage first."
            )

        # Check if outputs already exist
        if (
            self.fm.skip_if_exists(self._candidates_idlist, force)
            and self.fm.skip_if_exists(self._summary_tsv, force)
        ):
            return

        self.fm.ensure_dirs()

        # Step 1: determine target domains
        if self._domain.target_domains:
            target_domains: set[str] = set(self._domain.target_domains)
            target_superfamilies: set[str] = set()
            log.info("Using configured target domains: %s", target_domains)
        else:
            target_domains, target_superfamilies = self.extract_target_domains(self._cdd_result)
            log.info(
                "Auto-extracted %d domains, %d superfamilies",
                len(target_domains),
                len(target_superfamilies),
            )

        # Step 2: get genome CDD result
        if self._domain.genome_cdd and Path(self._domain.genome_cdd).exists():
            genome_cdd_path = Path(self._domain.genome_cdd)
            log.info("Using configured genome CDD: %s", genome_cdd_path)
        else:
            log.info("Submitting genome proteins to NCBI CD-Search...")
            pep_db = self.fm.rep_pep
            genome_cdd_path = self.submit_genome_batch(pep_db)

        # Step 3: parse genome CDD and filter
        genome_rows = self.parse_cdd_result(genome_cdd_path)
        filtered_ids = self.filter_by_domains(genome_rows, target_domains, target_superfamilies)
        log.info("Filtered %d genes by domain", len(filtered_ids))

        # Step 4: write candidates idlist
        self._candidates_idlist.write_text("\n".join(filtered_ids) + "\n")
        log.info("Candidates written → %s", self._candidates_idlist)

        # Step 5: write summary TSV
        filtered_set = set(filtered_ids)
        summary_rows = [
            row for row in genome_rows if row.gene_id in filtered_set
        ]
        lines = ["gene_id\tdomain_accession\tdomain_name\tsuperfamily_accession\tevalue"]
        for row in summary_rows:
            lines.append(
                f"{row.gene_id}\t{row.domain_accession}\t{row.domain_name}"
                f"\t{row.superfamily_accession}\t{row.evalue}"
            )
        self._summary_tsv.write_text("\n".join(lines) + "\n")
        log.info("Summary written → %s", self._summary_tsv)

        log.info("Domain-filter stage complete.")
