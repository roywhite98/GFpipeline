"""Unit tests for core/sequence.py."""

from pathlib import Path

import pytest

from gfpipeline.core.sequence import (
    extract_fasta,
    parse_blast_idlist,
    parse_hmm_idlist,
    transcript_to_gene_id,
)


# ---------------------------------------------------------------------------
# transcript_to_gene_id
# ---------------------------------------------------------------------------

class TestTranscriptToGeneId:
    def test_basic_conversion(self):
        assert transcript_to_gene_id("Os01t0936800-01") == "Os01g0936800"

    def test_different_isoform_suffix(self):
        assert transcript_to_gene_id("Os12t0123456-02") == "Os12g0123456"

    def test_same_gene_different_isoforms(self):
        g1 = transcript_to_gene_id("Os01t0936800-01")
        g2 = transcript_to_gene_id("Os01t0936800-02")
        assert g1 == g2

    def test_removes_isoform_suffix(self):
        result = transcript_to_gene_id("Os05t0400900-01")
        assert "-" not in result

    def test_replaces_t_with_g(self):
        result = transcript_to_gene_id("Os03t0285800-01")
        assert "g" in result
        assert result == "Os03g0285800"


# ---------------------------------------------------------------------------
# parse_hmm_idlist
# ---------------------------------------------------------------------------

HMM_OUT_CONTENT = """\
# hmmsearch :: search profile(s) against a sequence database
# HMMER 3.3.2 (Nov 2020); http://hmmer.org/
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# query HMM file:                  test.hmm
# target sequence database:        pep.fa
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

Query:       test  [M=200]
Scores for complete sequences (score includes all domains):
   --- full sequence ---   --- best 1 domain ---    -#dom-
    E-value  score  bias    E-value  score  bias    exp  N  Sequence       Description
    ------- ------ -----    ------- ------ -----   ---- --  --------       -----------
    1.2e-50  170.3   0.1    1.5e-50  170.0   0.1    1.0  1  Os01g0936800   gene1
    3.4e-10   38.2   0.0    4.0e-10   37.9   0.0    1.0  1  Os02g0123456   gene2
    0.0015    12.1   0.0    0.0020    11.8   0.0    1.0  1  Os03g0999999   gene3 (above threshold)
    0.5       5.0    0.0    0.6       4.8    0.0    1.0  1  Os04g0111111   gene4 (way above)

Domain annotation for each sequence (and alignments):
"""


class TestParseHmmIdlist:
    def test_extracts_ids_below_threshold(self, tmp_path):
        f = tmp_path / "hmm.out"
        f.write_text(HMM_OUT_CONTENT)
        ids = parse_hmm_idlist(f, evalue_threshold=100)
        assert "Os01g0936800" in ids
        assert "Os02g0123456" in ids
        assert "Os03g0999999" in ids

    def test_excludes_ids_above_threshold(self, tmp_path):
        f = tmp_path / "hmm.out"
        f.write_text(HMM_OUT_CONTENT)
        ids = parse_hmm_idlist(f, evalue_threshold=0.001)
        assert "Os01g0936800" in ids
        assert "Os02g0123456" in ids
        assert "Os03g0999999" not in ids
        assert "Os04g0111111" not in ids

    def test_returns_sorted_list(self, tmp_path):
        f = tmp_path / "hmm.out"
        f.write_text(HMM_OUT_CONTENT)
        ids = parse_hmm_idlist(f, evalue_threshold=100)
        assert ids == sorted(ids)

    def test_returns_deduplicated_list(self, tmp_path):
        f = tmp_path / "hmm.out"
        f.write_text(HMM_OUT_CONTENT)
        ids = parse_hmm_idlist(f, evalue_threshold=100)
        assert len(ids) == len(set(ids))

    def test_empty_result_when_threshold_zero(self, tmp_path):
        f = tmp_path / "hmm.out"
        f.write_text(HMM_OUT_CONTENT)
        ids = parse_hmm_idlist(f, evalue_threshold=0.0)
        assert ids == []


# ---------------------------------------------------------------------------
# parse_blast_idlist
# ---------------------------------------------------------------------------

BLAST_OUT_CONTENT = """\
Os01g0936800\tOs02g0111111\t95.0\t200\t10\t0\t1\t200\t1\t200\t1e-100\t300
Os01g0936800\tOs03g0222222\t80.0\t180\t20\t0\t1\t180\t1\t180\t1e-80\t250
Os05g0333333\tOs02g0111111\t70.0\t150\t30\t0\t1\t150\t1\t150\t1e-60\t200
Os05g0333333\tOs02g0111111\t65.0\t140\t35\t0\t1\t140\t1\t140\t1e-55\t190
"""


class TestParseBlastIdlist:
    def test_extracts_subject_ids(self, tmp_path):
        f = tmp_path / "blast.out"
        f.write_text(BLAST_OUT_CONTENT)
        ids = parse_blast_idlist(f)
        assert "Os02g0111111" in ids
        assert "Os03g0222222" in ids

    def test_deduplicates_ids(self, tmp_path):
        f = tmp_path / "blast.out"
        f.write_text(BLAST_OUT_CONTENT)
        ids = parse_blast_idlist(f)
        assert len(ids) == len(set(ids))
        # Os02g0111111 appears twice in input but once in output
        assert ids.count("Os02g0111111") == 1

    def test_returns_sorted_list(self, tmp_path):
        f = tmp_path / "blast.out"
        f.write_text(BLAST_OUT_CONTENT)
        ids = parse_blast_idlist(f)
        assert ids == sorted(ids)

    def test_skips_comment_lines(self, tmp_path):
        content = "# comment line\n" + BLAST_OUT_CONTENT
        f = tmp_path / "blast.out"
        f.write_text(content)
        ids = parse_blast_idlist(f)
        assert "# comment line" not in ids

    def test_empty_file_returns_empty_list(self, tmp_path):
        f = tmp_path / "blast.out"
        f.write_text("")
        assert parse_blast_idlist(f) == []


# ---------------------------------------------------------------------------
# extract_fasta
# ---------------------------------------------------------------------------

FASTA_CONTENT = """\
>Os01g0936800 gene1 description
MAAAKLLVVVGGGGG
>Os02g0123456 gene2 description
MTTTTLLLLLAAAAA
>Os03g0285800 gene3 description
MCCCCWWWWWDDDDD
"""


class TestExtractFasta:
    def test_extracts_requested_ids(self, tmp_path):
        db = tmp_path / "pep.fa"
        db.write_text(FASTA_CONTENT)
        out = tmp_path / "out.fa"
        count = extract_fasta(db, ["Os01g0936800", "Os03g0285800"], out)
        assert count == 2
        content = out.read_text()
        assert ">Os01g0936800" in content
        assert ">Os03g0285800" in content
        assert ">Os02g0123456" not in content

    def test_returns_count_of_extracted(self, tmp_path):
        db = tmp_path / "pep.fa"
        db.write_text(FASTA_CONTENT)
        out = tmp_path / "out.fa"
        count = extract_fasta(db, ["Os01g0936800"], out)
        assert count == 1

    def test_missing_ids_are_skipped(self, tmp_path):
        db = tmp_path / "pep.fa"
        db.write_text(FASTA_CONTENT)
        out = tmp_path / "out.fa"
        count = extract_fasta(db, ["Os01g0936800", "NONEXISTENT"], out)
        assert count == 1

    def test_empty_id_list_produces_empty_output(self, tmp_path):
        db = tmp_path / "pep.fa"
        db.write_text(FASTA_CONTENT)
        out = tmp_path / "out.fa"
        count = extract_fasta(db, [], out)
        assert count == 0
