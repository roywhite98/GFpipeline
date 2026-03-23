"""Domain stage: submit candidate proteins to NCBI Batch CD-Search."""

from __future__ import annotations

import logging
import time
from pathlib import Path

import requests
from Bio import SeqIO

from gfpipeline.config.schema import PipelineConfig
from gfpipeline.core.exceptions import ApiError, StageInputError
from gfpipeline.core.file_manager import FileManager

log = logging.getLogger(__name__)

_BWRPSB_URL = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"
_BATCH_SIZE = 4000
_POLL_INTERVAL = 5  # seconds
_BATCH_DELAY = 2    # seconds between batch submissions


class DomainStage:
    """Conserved domain analysis via NCBI Batch CD-Search.

    Steps:
        1. Check candidates.pep.fa exists
        2. Submit sequences in batches (≤4000 per batch)
        3. Poll until complete
        4. Download and save results to domain.cdd.txt
    """

    def __init__(self, config: PipelineConfig, fm: FileManager) -> None:
        self.config = config
        self.fm = fm
        self._domain = config.domain

    # ------------------------------------------------------------------
    # Path helpers
    # ------------------------------------------------------------------

    @property
    def _pep_fa(self) -> Path:
        return self.fm.result("identify", "candidates.pep", "fa")

    @property
    def _cdd_output(self) -> Path:
        return self.fm.result("domain", "cdd", "txt")

    # ------------------------------------------------------------------
    # API methods
    # ------------------------------------------------------------------

    def submit_batch(self, fasta_path: Path) -> str:
        """Submit sequences from fasta_path to NCBI Batch CD-Search, return RID."""
        fasta_content = fasta_path.read_text()
        data = {
            "useid1": "true",
            "maxhit": str(self._domain.maxhit),
            "filter": "true",
            "db": "cdd",
            "evalue": str(self._domain.evalue),
            "tdata": "hits",
            "queries": fasta_content,
        }
        resp = requests.post(_BWRPSB_URL, data=data)
        if resp.status_code != 200:
            raise ApiError(resp.status_code, resp.reason)

        # Parse RID from response text
        rid = self._parse_rid(resp.text)
        log.info("Submitted batch, RID=%s", rid)
        return rid

    def _parse_rid(self, text: str) -> str:
        """Extract CDSID/RID from API response text."""
        for line in text.splitlines():
            line = line.strip()
            if line.startswith("#cdsid"):
                # Format: #cdsid    QM3-qcdsearch-XXXXXXXX
                parts = line.split()
                if len(parts) >= 2:
                    return parts[1]
        raise ApiError(reason=f"Could not parse RID from response: {text[:200]}")

    def poll_status(self, rid: str) -> None:
        """Poll every 5 seconds until status code is 0; raise ApiError on error codes."""
        while True:
            print("等待结果中...")
            resp = requests.post(_BWRPSB_URL, data={"tdata": "hits", "cdsid": rid})
            if resp.status_code != 200:
                raise ApiError(resp.status_code, resp.reason)

            status = self._parse_status(resp.text)
            if status == 0:
                log.info("CD-Search job %s complete.", rid)
                return
            elif status == 3:
                time.sleep(_POLL_INTERVAL)
            elif status in (1, 2, 4, 5):
                raise ApiError(reason=f"CD-Search job {rid} failed with status code {status}")
            else:
                # Unknown status, keep polling
                time.sleep(_POLL_INTERVAL)

    def _parse_status(self, text: str) -> int:
        """Extract status code from API response text."""
        for line in text.splitlines():
            line = line.strip()
            if line.startswith("#status"):
                parts = line.split()
                if len(parts) >= 2:
                    return int(parts[1])
        # If no status line found, assume still in progress
        return 3

    def download_result(self, rid: str, output: Path) -> None:
        """Download CD-Search result and save to output path."""
        resp = requests.post(_BWRPSB_URL, data={"tdata": "hits", "cdsid": rid})
        if resp.status_code != 200:
            raise ApiError(resp.status_code, resp.reason)
        output.write_text(resp.text)
        log.info("Downloaded result to %s", output)

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------

    def run(self, force: bool = False) -> None:
        """Execute the full domain stage."""
        if not self._pep_fa.exists():
            raise StageInputError(
                f"Input file '{self._pep_fa}' does not exist. "
                "Please run the identify stage first."
            )

        if self.fm.skip_if_exists(self._cdd_output, force):
            return

        self.fm.ensure_dirs()

        # Read all sequences and split into batches
        records = list(SeqIO.parse(str(self._pep_fa), "fasta"))
        log.info("Total sequences: %d", len(records))

        all_results: list[str] = []
        header_written = False

        for i in range(0, len(records), _BATCH_SIZE):
            batch = records[i: i + _BATCH_SIZE]
            batch_num = i // _BATCH_SIZE + 1
            log.info("Submitting batch %d (%d sequences)", batch_num, len(batch))

            # Write batch to a temp fasta string
            import io
            buf = io.StringIO()
            SeqIO.write(batch, buf, "fasta")
            batch_fasta = buf.getvalue()

            # Submit
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

            # Poll
            self.poll_status(rid)

            # Download
            resp = requests.post(_BWRPSB_URL, data={"tdata": "hits", "cdsid": rid})
            if resp.status_code != 200:
                raise ApiError(resp.status_code, resp.reason)

            # Collect result lines (skip comment/header lines after first batch)
            for line in resp.text.splitlines():
                if line.startswith("#"):
                    if not header_written:
                        all_results.append(line)
                else:
                    all_results.append(line)
            header_written = True

            if i + _BATCH_SIZE < len(records):
                time.sleep(_BATCH_DELAY)

        self._cdd_output.write_text("\n".join(all_results) + "\n")
        log.info("Domain stage complete → %s", self._cdd_output)
