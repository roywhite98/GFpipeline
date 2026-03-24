"""Domain-filter stage: filter genome proteins by conserved domain."""

from __future__ import annotations

import io
import logging
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import requests
import click
from Bio import SeqIO

from gfpipeline.config.schema import PipelineConfig
from gfpipeline.core.exceptions import ApiError, StageInputError
from gfpipeline.core.file_manager import FileManager

log = logging.getLogger(__name__)

_BWRPSB_URL = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"
_BATCH_SIZE = 4000
_POLL_INTERVAL = 5
_BATCH_DELAY = 2


def _is_placeholder(value: str) -> bool:
    """Return True if value is an empty/placeholder string (e.g. '-', ' - ', '')."""
    return not value or value.strip() in ("", "-")


# ---------------------------------------------------------------------------
# Boolean query parser
# ---------------------------------------------------------------------------
# Grammar (AND/OR/NOT are case-insensitive keywords):
#   expr   := term (OR term)*
#   term   := factor (AND factor)*
#   factor := NOT factor | '(' expr ')' | TOKEN
#
# Returns a callable: (gene_domains: set[str]) -> bool
#
# Example queries:
#   "pfam00161 OR cl08249"
#   "(pfam00161 OR cl08249) AND pfam20241"
#   "pfam00161 AND NOT pfam20241"

_TOKEN_RE = re.compile(r'\(|\)|[^\s()]+')


def _tokenize(query: str) -> list[str]:
    return _TOKEN_RE.findall(query)


class _Parser:
    """Recursive-descent boolean expression parser."""

    def __init__(self, tokens: list[str]) -> None:
        self._tokens = tokens
        self._pos = 0

    def _peek(self) -> str | None:
        return self._tokens[self._pos] if self._pos < len(self._tokens) else None

    def _consume(self) -> str:
        tok = self._tokens[self._pos]
        self._pos += 1
        return tok

    def parse(self) -> Callable[[set[str]], bool]:
        node = self._expr()
        if self._peek() is not None:
            raise ValueError(f"Unexpected token '{self._peek()}' in domain query")
        return node

    def _expr(self) -> Callable[[set[str]], bool]:
        left = self._term()
        while self._peek() and self._peek().upper() == "OR":
            self._consume()
            right = self._term()
            left = (lambda l, r: lambda d: l(d) or r(d))(left, right)
        return left

    def _term(self) -> Callable[[set[str]], bool]:
        left = self._factor()
        while self._peek() and self._peek().upper() == "AND":
            self._consume()
            right = self._factor()
            left = (lambda l, r: lambda d: l(d) and r(d))(left, right)
        return left

    def _factor(self) -> Callable[[set[str]], bool]:
        tok = self._peek()
        if tok is None:
            raise ValueError("Unexpected end of domain query")
        if tok.upper() == "NOT":
            self._consume()
            inner = self._factor()
            return lambda d, f=inner: not f(d)
        if tok == "(":
            self._consume()
            node = self._expr()
            if self._peek() != ")":
                raise ValueError("Missing closing ')' in domain query")
            self._consume()
            return node
        self._consume()
        return lambda d, t=tok: t in d


def parse_domain_query(query: str) -> Callable[[set[str]], bool]:
    """Parse a boolean domain query string into an evaluator function.

    Supports AND, OR, NOT (case-insensitive) and parentheses for grouping.

    Examples:
        "pfam00161"                              -> matches genes with pfam00161
        "pfam00161 OR cl08249"                   -> matches either
        "(pfam00161 OR cl08249) AND pfam20241"   -> must have pfam20241 plus one of the others
        "pfam00161 AND NOT pfam20241"            -> has pfam00161 but not pfam20241

    Returns a callable: (gene_domains: set[str]) -> bool
    Raises ValueError on syntax errors.
    """
    tokens = _tokenize(query.strip())
    if not tokens:
        raise ValueError("Empty domain query")
    return _Parser(tokens).parse()


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class DomainSummaryRow:
    gene_id: str
    domain_accession: str
    domain_name: str
    superfamily_accession: str
    evalue: float


# ---------------------------------------------------------------------------
# Stage
# ---------------------------------------------------------------------------

class DomainFilterStage:
    """Filter genome proteins by conserved domain using a boolean query.

    Config field ``domain.target_domains`` accepts a boolean query string, e.g.:
        "pfam00161 OR cl08249"
        "(pfam00161 OR cl08249) AND pfam20241"

    If not set, the intersection of domains shared by all ref sequences is used
    as the auto-extracted target (AND logic across all ref genes).
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
    def _cdd_dir(self) -> Path:
        """data/cdd/ directory for genome-level CDD files."""
        return self.fm.data_dir / "cdd"

    @property
    def _genome_cdd(self) -> Path:
        genome_name = self.config.genome_db.genome_name or self.config.project_name
        return self._cdd_dir / f"{genome_name}.genome.cdd.txt"

    @property
    def _candidates_idlist(self) -> Path:
        return self.fm.result("domain-filter", "candidates", "idlist")

    @property
    def _summary_tsv(self) -> Path:
        return self.fm.result("domain-filter", "summary", "tsv")

    # ------------------------------------------------------------------
    # Core methods
    # ------------------------------------------------------------------

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
            if parts[0].strip() == "Query":
                continue
            query = parts[0].strip()
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

    def filter_by_query(
        self,
        rows: list[DomainSummaryRow],
        evaluator: Callable[[set[str]], bool],
    ) -> list[str]:
        """Return sorted list of gene IDs whose domain set satisfies the boolean evaluator."""
        gene_domains: dict[str, set[str]] = {}
        for row in rows:
            if _is_placeholder(row.domain_accession):
                continue
            gene_domains.setdefault(row.gene_id, set()).add(row.domain_accession)

        return sorted(
            gene_id for gene_id, domains in gene_domains.items()
            if evaluator(domains)
        )

    def auto_extract_query(self, cdd_result: Path) -> Callable[[set[str]], bool] | None:
        """Build an AND evaluator from the intersection of domains in all ref sequences.

        Returns None if the intersection is empty (no shared domain found).
        """
        rows = self.parse_cdd_result(cdd_result)
        gene_domains: dict[str, set[str]] = {}
        for row in rows:
            if _is_placeholder(row.domain_accession):
                continue
            gene_domains.setdefault(row.gene_id, set()).add(row.domain_accession)

        if not gene_domains:
            log.warning("auto_extract_query: no valid domain hits found in ref CDD result")
            return None

        intersection = set.intersection(*gene_domains.values())
        if not intersection:
            log.warning(
                "auto_extract_query: domain intersection is empty "
                "(no domain shared by all %d ref sequences) — "
                "please set domain.target_domains in config.yaml",
                len(gene_domains),
            )
            return None

        log.info("Auto-extracted domain intersection: %s (from %d ref sequences)", intersection, len(gene_domains))
        # AND of all shared domains
        def _evaluator(domains: set[str], required: set[str] = intersection) -> bool:
            return required.issubset(domains)
        return _evaluator

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
    # Local rpsblast
    # ------------------------------------------------------------------

    def run_local_rpsblast(self, pep_db: Path) -> Path:
        """Run rpsblast locally against a CDD database, convert output to bwrpsb format."""
        import subprocess
        cdd_db = self._domain.cdd_db
        rpsblast = self.config.tools.rpsblast
        raw_out = self._genome_cdd.with_suffix(".rpsblast.tsv")

        log.info("Running local rpsblast: %s -db %s -query %s", rpsblast, cdd_db, pep_db)
        cmd = [
            rpsblast,
            "-query", str(pep_db),
            "-db", cdd_db,
            "-out", str(raw_out),
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle",
            "-evalue", str(self._domain.evalue),
            "-num_threads", "10",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"rpsblast failed:\n{result.stderr}")
        log.info("rpsblast done → %s", raw_out)

        lines = ["#Query\tHit_type\tPSSM_ID\tFrom\tTo\tE-Value\tBitscore\tAccession\tShort_name\tIncomplete\tSuperfamily"]
        for line in raw_out.read_text().splitlines():
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 13:
                continue
            qseqid   = parts[0]
            sseqid   = parts[1]
            evalue   = parts[10]
            bitscore = parts[11]
            stitle   = parts[12]
            title_parts = stitle.split(",", 2)
            accession  = title_parts[0].strip() if len(title_parts) > 0 else sseqid
            short_name = title_parts[1].strip() if len(title_parts) > 1 else ""
            pssm_id = sseqid.split("|")[-1] if "|" in sseqid else sseqid
            lines.append(
                f"{qseqid}\tspecific\t{pssm_id}\t-\t-\t{evalue}\t{bitscore}\t{accession}\t{short_name}\t-\t-"
            )

        self._genome_cdd.write_text("\n".join(lines) + "\n")
        log.info("Converted CDD result saved → %s", self._genome_cdd)
        return self._genome_cdd

    # ------------------------------------------------------------------
    # API helpers
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

        if (
            self.fm.skip_if_exists(self._candidates_idlist, force)
            and self.fm.skip_if_exists(self._summary_tsv, force)
        ):
            return

        self.fm.ensure_dirs()
        self._cdd_dir.mkdir(parents=True, exist_ok=True)

        # Step 1: build evaluator
        if self._domain.target_domains:
            query_str = self._domain.target_domains
            evaluator = parse_domain_query(query_str)
            log.info("Using domain query: %s", query_str)
        else:
            evaluator = self.auto_extract_query(self._cdd_result)
            if evaluator is None:
                raise StageInputError(
                    "Could not determine target domains automatically. "
                    "Please set domain.target_domains in config.yaml, "
                    "e.g.: target_domains: \"pfam00161 OR cl08249\""
                )

        # Step 2: get genome CDD result
        if self._domain.genome_cdd and Path(self._domain.genome_cdd).exists():
            genome_cdd_path = Path(self._domain.genome_cdd)
            log.info("Using configured genome CDD: %s", genome_cdd_path)
        elif self._domain.cdd_db:
            log.info("Running local rpsblast with CDD db: %s", self._domain.cdd_db)
            click.echo("正在运行全基因组 rpsblast，序列较多时耗时较长，请耐心等待...")
            genome_cdd_path = self.run_local_rpsblast(self.fm.rep_pep)
        else:
            log.info("Submitting genome proteins to NCBI CD-Search...")
            click.echo("正在向 NCBI CD-Search 提交全基因组蛋白，耗时较长（通常数分钟至数十分钟），请耐心等待...")
            genome_cdd_path = self.submit_genome_batch(self.fm.rep_pep)

        # Step 3: filter
        genome_rows = self.parse_cdd_result(genome_cdd_path)
        filtered_ids = self.filter_by_query(genome_rows, evaluator)
        log.info("Filtered %d genes by domain query", len(filtered_ids))

        # Step 4: write candidates idlist
        self._candidates_idlist.write_text("\n".join(filtered_ids) + "\n")
        log.info("Candidates written → %s", self._candidates_idlist)

        # Step 5: write summary TSV
        filtered_set = set(filtered_ids)
        summary_rows = [row for row in genome_rows if row.gene_id in filtered_set]
        lines = ["gene_id\tdomain_accession\tdomain_name\tsuperfamily_accession\tevalue"]
        for row in summary_rows:
            lines.append(
                f"{row.gene_id}\t{row.domain_accession}\t{row.domain_name}"
                f"\t{row.superfamily_accession}\t{row.evalue}"
            )
        self._summary_tsv.write_text("\n".join(lines) + "\n")
        log.info("Summary written → %s", self._summary_tsv)

        log.info("Domain-filter stage complete.")
