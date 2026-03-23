"""Representative transcript index builder for genome-db stage.

Selects one representative transcript per gene, extracts CDS sequences,
translates to protein, and writes output files.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Optional

from gfpipeline.config.schema import PipelineConfig
from gfpipeline.core.exceptions import PipelineError, StageInputError
from gfpipeline.core.file_manager import FileManager
from gfpipeline.core.logger import get_logger

log = get_logger(__name__)


def _parse_attributes(attr_str: str) -> dict[str, str]:
    """Parse GFF3 attribute column into a dict."""
    attrs: dict[str, str] = {}
    for part in attr_str.strip().split(";"):
        part = part.strip()
        if "=" in part:
            k, _, v = part.partition("=")
            attrs[k.strip()] = v.strip()
    return attrs


class RepIndexBuilder:
    """Select representative transcripts and build rep sequence indices.

    Args:
        config: Pipeline configuration.
        fm: File manager.
    """

    def __init__(self, config: PipelineConfig, fm: FileManager) -> None:
        self.config = config
        self.fm = fm

    # ------------------------------------------------------------------
    # Public methods
    # ------------------------------------------------------------------

    def select_representative(
        self,
        gene_id: str,
        transcripts: list[str],
        cds_lengths: dict[str, int],
    ) -> Optional[str]:
        """Select the representative transcript for a gene.

        Strategy 'longest_cds' (default): pick transcript with longest CDS;
        ties broken by lexicographic order (smallest ID wins).

        Strategy 'longest_mrna': pick transcript with longest mRNA span;
        ties broken by lexicographic order.

        Args:
            gene_id: Gene identifier (used for logging).
            transcripts: List of transcript IDs for this gene.
            cds_lengths: Mapping of transcript_id -> CDS length (or mRNA length).

        Returns:
            Selected transcript ID, or None if no valid transcript.
        """
        if not transcripts:
            log.warning("%s: no transcripts, skipping", gene_id)
            return None

        valid = {t: cds_lengths.get(t, 0) for t in transcripts if cds_lengths.get(t, 0) > 0}
        if not valid:
            log.warning("%s: all transcripts have no CDS, skipping", gene_id)
            return None

        max_len = max(valid.values())
        candidates = sorted(t for t, length in valid.items() if length == max_len)
        return candidates[0]  # lexicographically smallest

    def run(self, force: bool = False) -> None:
        """Build representative transcript index.

        Steps:
        1. Parse GFF3 to get gene->transcript mapping and CDS/mRNA lengths.
        2. Select representative transcript per gene.
        3. Extract CDS sequences -> rep.cds.fa.
        4. Translate to protein -> rep.pep.fa.
        5. Write rep.idlist and gene2rep.tsv.

        Args:
            force: If True, rebuild even if outputs already exist.

        Raises:
            StageInputError: If required input files are missing.
        """
        db_cfg = self.config.databases
        cfg = self.config.genome_db

        gff3_path = Path(db_cfg.gff3)
        genome_path = Path(db_cfg.genome)

        for label, p in [("GFF3", gff3_path), ("genome", genome_path)]:
            if not p.exists():
                raise StageInputError(f"rep-index: {label} file not found: {p}")

        genome_name = cfg.genome_name or self.config.project_name
        rep_dir = Path(cfg.rep_index_dir)
        rep_dir.mkdir(parents=True, exist_ok=True)

        cds_out = rep_dir / f"{genome_name}.rep.cds.fa"
        pep_out = rep_dir / f"{genome_name}.rep.pep.fa"
        idlist_out = rep_dir / f"{genome_name}.rep.idlist"
        g2rep_out = rep_dir / f"{genome_name}.gene2rep.tsv"

        if not force and cds_out.exists() and g2rep_out.exists():
            log.info("Skipping rep-index (outputs already exist). Use --force to re-run.")
            return

        # Parse GFF3
        gene2transcripts, cds_lengths, mrna_lengths = self._parse_gff3(gff3_path)

        # Select lengths based on strategy
        if cfg.rep_selection == "longest_mrna":
            lengths_to_use = mrna_lengths
        else:
            lengths_to_use = cds_lengths

        # Select representatives
        gene2rep: dict[str, str] = {}
        for gene_id, transcripts in sorted(gene2transcripts.items()):
            rep = self.select_representative(gene_id, transcripts, lengths_to_use)
            if rep is not None:
                gene2rep[gene_id] = rep

        log.info("Selected %d representative transcripts", len(gene2rep))

        # Extract CDS sequences from genome + GFF3
        rep_ids = list(gene2rep.values())
        self._extract_cds_from_genome(genome_path, gff3_path, rep_ids, cds_out)

        # Translate to protein
        self._translate_cds(cds_out, pep_out)

        # Write idlist
        with open(idlist_out, "w") as fh:
            for t_id in sorted(rep_ids):
                fh.write(t_id + "\n")
        log.info("Wrote rep.idlist: %s", idlist_out)

        # Write gene2rep.tsv
        with open(g2rep_out, "w") as fh:
            fh.write("gene_id\trep_transcript_id\tcds_length\ttranscript_count\n")
            for gene_id in sorted(gene2rep.keys()):
                rep_t = gene2rep[gene_id]
                cds_len = cds_lengths.get(rep_t, 0)
                t_count = len(gene2transcripts.get(gene_id, []))
                fh.write(f"{gene_id}\t{rep_t}\t{cds_len}\t{t_count}\n")
        log.info("Wrote gene2rep.tsv: %s", g2rep_out)

        log.info("rep-index complete. Outputs in %s", rep_dir)

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _parse_gff3(
        self, gff3_path: Path
    ) -> tuple[dict[str, list[str]], dict[str, int], dict[str, int]]:
        """Parse GFF3 to extract gene->transcript mapping and lengths.

        Returns:
            Tuple of:
            - gene2transcripts: gene_id -> [transcript_id, ...]
            - cds_lengths: transcript_id -> total CDS length
            - mrna_lengths: transcript_id -> mRNA span length (end - start + 1)
        """
        gene2transcripts: dict[str, list[str]] = defaultdict(list)
        cds_lengths: dict[str, int] = defaultdict(int)
        mrna_lengths: dict[str, int] = {}
        transcript_to_gene: dict[str, str] = {}

        try:
            with open(gff3_path) as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.split("\t")
                    if len(parts) < 9:
                        continue
                    chrom, _, feat_type, start_s, end_s, _, strand, _, attr_s = parts[:9]
                    try:
                        start = int(start_s)
                        end = int(end_s)
                    except ValueError:
                        continue
                    attrs = _parse_attributes(attr_s)

                    if feat_type in ("mRNA", "transcript"):
                        t_id = attrs.get("ID", "")
                        parent = attrs.get("Parent", "")
                        if t_id and parent:
                            gene2transcripts[parent].append(t_id)
                            transcript_to_gene[t_id] = parent
                            mrna_lengths[t_id] = end - start + 1

                    elif feat_type == "CDS":
                        parent = attrs.get("Parent", "")
                        if parent:
                            cds_lengths[parent] += end - start + 1

        except OSError as exc:
            raise PipelineError(f"rep-index: cannot read GFF3 {gff3_path}: {exc}") from exc

        return dict(gene2transcripts), dict(cds_lengths), mrna_lengths

    def _extract_sequences(
        self, fasta_path: Path, ids: list[str], output_path: Path
    ) -> None:
        """Extract sequences for given IDs from a FASTA file using BioPython."""
        try:
            from Bio import SeqIO
        except ImportError:
            raise PipelineError("BioPython is required for sequence extraction.")

        id_set = set(ids)
        written = 0

        with open(output_path, "w") as fh:
            for record in SeqIO.parse(str(fasta_path), "fasta"):
                if record.id in id_set:
                    fh.write(f">{record.id}\n{str(record.seq)}\n")
                    written += 1

        log.info("Extracted %d sequences to %s", written, output_path)

    def _extract_cds_from_genome(
        self,
        genome_path: Path,
        gff3_path: Path,
        transcript_ids: list[str],
        output_path: Path,
    ) -> None:
        """Extract and splice CDS sequences for given transcript IDs from genome + GFF3.

        Parses GFF3 for CDS intervals, fetches genomic sequence, splices and
        reverse-complements as needed, then writes to output_path.

        Args:
            genome_path: Path to genome FASTA.
            gff3_path: Path to GFF3 annotation file.
            transcript_ids: List of transcript IDs to extract.
            output_path: Output FASTA path.
        """
        try:
            from Bio import SeqIO
            from Bio.Seq import Seq
        except ImportError:
            raise PipelineError("BioPython is required for CDS extraction.")

        target_ids = set(transcript_ids)

        # Parse GFF3: collect CDS intervals per transcript
        # cds_map: transcript_id -> [(chrom, start, end, strand), ...]  (1-based coords)
        cds_map: dict[str, list[tuple[str, int, int, str]]] = defaultdict(list)
        with open(gff3_path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split("\t")
                if len(parts) < 9:
                    continue
                chrom, _, feat_type, start_s, end_s, _, strand, _, attr_s = parts[:9]
                if feat_type != "CDS":
                    continue
                try:
                    start = int(start_s)
                    end = int(end_s)
                except ValueError:
                    continue
                attrs = _parse_attributes(attr_s)
                parent = attrs.get("Parent", "")
                if parent in target_ids:
                    cds_map[parent].append((chrom, start, end, strand))

        # Load genome into memory (indexed by chrom id)
        genome_index = SeqIO.to_dict(SeqIO.parse(str(genome_path), "fasta"))

        written = 0
        with open(output_path, "w") as fh:
            for t_id in transcript_ids:
                intervals = cds_map.get(t_id)
                if not intervals:
                    log.warning("rep-index: no CDS intervals found for %s, skipping", t_id)
                    continue
                chrom = intervals[0][0]
                strand = intervals[0][3]
                if chrom not in genome_index:
                    log.warning("rep-index: chromosome '%s' not in genome, skipping %s", chrom, t_id)
                    continue
                chrom_seq = genome_index[chrom].seq
                # Sort CDS intervals by start position
                sorted_intervals = sorted(intervals, key=lambda x: x[1])
                cds_seq = Seq("".join(
                    str(chrom_seq[s - 1:e])  # convert 1-based to 0-based
                    for _, s, e, _ in sorted_intervals
                ))
                if strand == "-":
                    cds_seq = cds_seq.reverse_complement()
                fh.write(f">{t_id}\n{str(cds_seq)}\n")
                written += 1

        log.info("Extracted %d CDS sequences from genome to %s", written, output_path)

    def _translate_cds(self, cds_path: Path, pep_path: Path) -> None:
        """Translate CDS sequences to protein using BioPython."""
        try:
            from Bio import SeqIO
            from Bio.Seq import Seq
        except ImportError:
            raise PipelineError("BioPython is required for CDS translation.")

        written = 0
        with open(pep_path, "w") as fh:
            for record in SeqIO.parse(str(cds_path), "fasta"):
                seq = str(record.seq)
                # Trim to multiple of 3
                trim = len(seq) - (len(seq) % 3)
                seq = seq[:trim]
                try:
                    pep = str(Seq(seq).translate(to_stop=True))
                except Exception as exc:
                    log.warning("Translation failed for %s: %s", record.id, exc)
                    pep = ""
                if pep:
                    fh.write(f">{record.id}\n{pep}\n")
                    written += 1

        log.info("Translated %d sequences to %s", written, pep_path)
