"""Tests for DomainFilterStage."""

from __future__ import annotations

import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from gfpipeline.config.schema import DatabasesConfig, DomainConfig, GenomeDbConfig, PipelineConfig
from gfpipeline.core.exceptions import ApiError, StageInputError
from gfpipeline.core.file_manager import FileManager
from gfpipeline.stages.domain_filter import DomainFilterStage, DomainSummaryRow


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_config(
    tmp_path: Path,
    domain_cfg: DomainConfig | None = None,
) -> PipelineConfig:
    return PipelineConfig(
        project_name="TEST",
        data_dir=str(tmp_path / "data"),
        result_dir=str(tmp_path / "results"),
        databases=DatabasesConfig(
            ref_fa=str(tmp_path / "ref.fa"),
            genome=str(tmp_path / "genome.fa"),
            gff3=str(tmp_path / "anno.gff3"),
        ),
        domain=domain_cfg or DomainConfig(),
        genome_db=GenomeDbConfig(rep_index_dir=str(tmp_path / "rep_index")),
    )


def _make_stage(tmp_path: Path, domain_cfg: DomainConfig | None = None) -> DomainFilterStage:
    config = _make_config(tmp_path, domain_cfg)
    fm = FileManager(config)
    return DomainFilterStage(config, fm)


_SAMPLE_CDD = textwrap.dedent("""\
    #Batch CD-search tool
    #cdsid    QM3-qcdsearch-12345678-90AB
    Query\tHit_type\tPSSM_ID\tFrom\tTo\tE-Value\tBitscore\tAccession\tShort_name\tIncomplete\tSuperfamily
    gene1\tspecific\t12345\t1\t100\t1e-10\t200\tcd12345\tDomain_A\t-\tcl99999
    gene2\tspecific\t12346\t1\t80\t1e-5\t150\tcd99999\tDomain_B\t-\tcl11111
    gene3\tspecific\t12347\t5\t90\t0.001\t100\tcd00001\tDomain_C\t-\tcl99999
""")


# ---------------------------------------------------------------------------
# Tests: parse_cdd_result
# ---------------------------------------------------------------------------

def test_parse_cdd_result_basic(tmp_path):
    stage = _make_stage(tmp_path)
    cdd_file = tmp_path / "test.cdd.txt"
    cdd_file.write_text(_SAMPLE_CDD)

    rows = stage.parse_cdd_result(cdd_file)
    assert len(rows) == 3
    assert rows[0].gene_id == "gene1"
    assert rows[0].domain_accession == "cd12345"
    assert rows[0].domain_name == "Domain_A"
    assert rows[0].superfamily_accession == "cl99999"
    assert rows[0].evalue == pytest.approx(1e-10)


def test_parse_cdd_result_skips_comments(tmp_path):
    stage = _make_stage(tmp_path)
    cdd_file = tmp_path / "test.cdd.txt"
    cdd_file.write_text("# comment\n# another comment\n")
    rows = stage.parse_cdd_result(cdd_file)
    assert rows == []


def test_parse_cdd_result_missing_file_raises(tmp_path):
    stage = _make_stage(tmp_path)
    with pytest.raises(StageInputError):
        stage.parse_cdd_result(tmp_path / "nonexistent.txt")


def test_parse_cdd_result_skips_short_lines(tmp_path):
    stage = _make_stage(tmp_path)
    cdd_file = tmp_path / "test.cdd.txt"
    # Only 5 columns, should be skipped
    cdd_file.write_text("gene1\tspecific\t12345\t1\t100\n")
    rows = stage.parse_cdd_result(cdd_file)
    assert rows == []


# ---------------------------------------------------------------------------
# Tests: extract_target_domains
# ---------------------------------------------------------------------------

def test_extract_target_domains(tmp_path):
    stage = _make_stage(tmp_path)
    cdd_file = tmp_path / "test.cdd.txt"
    cdd_file.write_text(_SAMPLE_CDD)

    domains, superfamilies = stage.extract_target_domains(cdd_file)
    assert "cd12345" in domains
    assert "cd99999" in domains
    assert "cd00001" in domains
    assert "cl99999" in superfamilies
    assert "cl11111" in superfamilies


def test_extract_target_domains_empty_file(tmp_path):
    stage = _make_stage(tmp_path)
    cdd_file = tmp_path / "empty.cdd.txt"
    cdd_file.write_text("# only comments\n")

    domains, superfamilies = stage.extract_target_domains(cdd_file)
    assert domains == set()
    assert superfamilies == set()


# ---------------------------------------------------------------------------
# Tests: filter_by_domains
# ---------------------------------------------------------------------------

def _make_rows() -> list[DomainSummaryRow]:
    return [
        DomainSummaryRow("geneA", "cd12345", "Domain_A", "cl99999", 1e-10),
        DomainSummaryRow("geneB", "cd99999", "Domain_B", "cl11111", 1e-5),
        DomainSummaryRow("geneC", "cd00001", "Domain_C", "cl99999", 0.001),
        DomainSummaryRow("geneD", "cd77777", "Domain_D", "cl00000", 0.1),
    ]


def test_filter_by_domains_exact_match(tmp_path):
    stage = _make_stage(tmp_path)
    rows = _make_rows()
    result = stage.filter_by_domains(rows, {"cd12345"}, set())
    assert result == ["geneA"]


def test_filter_by_domains_superfamily_match(tmp_path):
    stage = _make_stage(tmp_path)
    rows = _make_rows()
    result = stage.filter_by_domains(rows, set(), {"cl99999"})
    # geneA and geneC both have cl99999 superfamily
    assert result == ["geneA", "geneC"]


def test_filter_by_domains_combined(tmp_path):
    stage = _make_stage(tmp_path)
    rows = _make_rows()
    result = stage.filter_by_domains(rows, {"cd99999"}, {"cl99999"})
    # geneA (superfamily cl99999), geneB (domain cd99999), geneC (superfamily cl99999)
    assert result == ["geneA", "geneB", "geneC"]


def test_filter_by_domains_no_match(tmp_path):
    stage = _make_stage(tmp_path)
    rows = _make_rows()
    result = stage.filter_by_domains(rows, {"cd_nonexistent"}, {"cl_nonexistent"})
    assert result == []


def test_filter_by_domains_returns_sorted(tmp_path):
    stage = _make_stage(tmp_path)
    rows = [
        DomainSummaryRow("geneZ", "cd12345", "Domain_A", "cl99999", 1e-10),
        DomainSummaryRow("geneA", "cd12345", "Domain_A", "cl99999", 1e-10),
        DomainSummaryRow("geneM", "cd12345", "Domain_A", "cl99999", 1e-10),
    ]
    result = stage.filter_by_domains(rows, {"cd12345"}, set())
    assert result == sorted(result)


def test_filter_by_domains_no_duplicates(tmp_path):
    stage = _make_stage(tmp_path)
    # Same gene appears multiple times with matching domain
    rows = [
        DomainSummaryRow("geneA", "cd12345", "Domain_A", "cl99999", 1e-10),
        DomainSummaryRow("geneA", "cd12345", "Domain_A", "cl99999", 1e-8),
    ]
    result = stage.filter_by_domains(rows, {"cd12345"}, set())
    assert result.count("geneA") == 1


def test_filter_by_domains_domain_takes_priority_over_superfamily(tmp_path):
    """Gene with matching domain should be included even if superfamily doesn't match."""
    stage = _make_stage(tmp_path)
    rows = [
        DomainSummaryRow("geneA", "cd12345", "Domain_A", "cl_other", 1e-10),
    ]
    result = stage.filter_by_domains(rows, {"cd12345"}, {"cl99999"})
    assert "geneA" in result


# ---------------------------------------------------------------------------
# Tests: run - input validation
# ---------------------------------------------------------------------------

def test_run_raises_if_cdd_missing(tmp_path):
    stage = _make_stage(tmp_path)
    with pytest.raises(StageInputError):
        stage.run()


def test_run_skips_if_outputs_exist(tmp_path):
    config = _make_config(tmp_path)
    fm = FileManager(config)
    fm.ensure_dirs()
    stage = DomainFilterStage(config, fm)

    # Create required input and outputs
    stage._cdd_result.write_text(_SAMPLE_CDD)
    stage._candidates_idlist.write_text("geneA\n")
    stage._summary_tsv.write_text("gene_id\tdomain_accession\n")

    with patch("gfpipeline.stages.domain_filter.requests.post") as mock_post:
        stage.run(force=False)
        mock_post.assert_not_called()


# ---------------------------------------------------------------------------
# Tests: run - with configured target_domains
# ---------------------------------------------------------------------------

def test_run_uses_configured_target_domains(tmp_path):
    """When target_domains is configured, skip auto-extraction."""
    domain_cfg = DomainConfig(target_domains=["cd12345"])
    config = _make_config(tmp_path, domain_cfg)
    fm = FileManager(config)
    fm.ensure_dirs()
    stage = DomainFilterStage(config, fm)

    # Create candidate CDD result
    stage._cdd_result.write_text(_SAMPLE_CDD)

    # Create genome CDD result (use same sample)
    genome_cdd = tmp_path / "genome.cdd.txt"
    genome_cdd.write_text(_SAMPLE_CDD)

    domain_cfg_with_genome = DomainConfig(
        target_domains=["cd12345"],
        genome_cdd=str(genome_cdd),
    )
    config2 = _make_config(tmp_path, domain_cfg_with_genome)
    fm2 = FileManager(config2)
    fm2.ensure_dirs()
    stage2 = DomainFilterStage(config2, fm2)
    stage2._cdd_result.write_text(_SAMPLE_CDD)

    with patch("gfpipeline.stages.domain_filter.requests.post") as mock_post:
        stage2.run()
        mock_post.assert_not_called()

    assert stage2._candidates_idlist.exists()
    ids = stage2._candidates_idlist.read_text().strip().splitlines()
    assert "gene1" in ids


# ---------------------------------------------------------------------------
# Tests: run - with genome_cdd configured
# ---------------------------------------------------------------------------

def test_run_uses_configured_genome_cdd(tmp_path):
    """When genome_cdd is configured and exists, skip API submission."""
    genome_cdd = tmp_path / "genome.cdd.txt"
    genome_cdd.write_text(_SAMPLE_CDD)

    domain_cfg = DomainConfig(genome_cdd=str(genome_cdd))
    config = _make_config(tmp_path, domain_cfg)
    fm = FileManager(config)
    fm.ensure_dirs()
    stage = DomainFilterStage(config, fm)
    stage._cdd_result.write_text(_SAMPLE_CDD)

    with patch("gfpipeline.stages.domain_filter.requests.post") as mock_post:
        stage.run()
        mock_post.assert_not_called()

    assert stage._candidates_idlist.exists()
    assert stage._summary_tsv.exists()


# ---------------------------------------------------------------------------
# Tests: run - full flow with mocked API
# ---------------------------------------------------------------------------

def _fake_submit_resp() -> MagicMock:
    return MagicMock(
        status_code=200,
        text="#cdsid    QM3-qcdsearch-12345678-90AB\n#status    3\n",
    )


def _fake_complete_resp() -> MagicMock:
    return MagicMock(status_code=200, text="#status    0\n")


def _fake_download_resp() -> MagicMock:
    return MagicMock(status_code=200, text=_SAMPLE_CDD)


def test_run_submits_genome_batch_when_no_genome_cdd(tmp_path):
    """When genome_cdd not configured, submit genome proteins via API."""
    config = _make_config(tmp_path)
    fm = FileManager(config)
    fm.ensure_dirs()
    stage = DomainFilterStage(config, fm)

    stage._cdd_result.write_text(_SAMPLE_CDD)
    # Create rep_pep for genome batch submission
    fm.rep_pep.parent.mkdir(parents=True, exist_ok=True)
    fm.rep_pep.write_text(">gene1\nMACDEFGH\n>gene2\nMKLMNOP\n")

    with patch(
        "gfpipeline.stages.domain_filter.requests.post",
        side_effect=[_fake_submit_resp(), _fake_complete_resp(), _fake_download_resp()],
    ):
        stage.run()

    assert stage._candidates_idlist.exists()
    assert stage._summary_tsv.exists()


# ---------------------------------------------------------------------------
# Tests: summary TSV format
# ---------------------------------------------------------------------------

def test_run_summary_tsv_has_header(tmp_path):
    genome_cdd = tmp_path / "genome.cdd.txt"
    genome_cdd.write_text(_SAMPLE_CDD)

    domain_cfg = DomainConfig(genome_cdd=str(genome_cdd))
    config = _make_config(tmp_path, domain_cfg)
    fm = FileManager(config)
    fm.ensure_dirs()
    stage = DomainFilterStage(config, fm)
    stage._cdd_result.write_text(_SAMPLE_CDD)

    stage.run()

    header = stage._summary_tsv.read_text().splitlines()[0]
    assert "gene_id" in header
    assert "domain_accession" in header
    assert "evalue" in header
