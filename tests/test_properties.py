"""Tests for PropertiesStage."""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from gfpipeline.config.schema import DatabasesConfig, PipelineConfig
from gfpipeline.core.exceptions import StageInputError
from gfpipeline.core.file_manager import FileManager
from gfpipeline.stages.properties import PropertiesStage, _AA_ORDER, _SCALAR_COLS


# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------

def _make_config(tmp_path: Path) -> PipelineConfig:
    return PipelineConfig(
        project_name="TEST",
        data_dir=str(tmp_path / "data"),
        result_dir=str(tmp_path / "results"),
        databases=DatabasesConfig(
            ref_fa=str(tmp_path / "ref.fa"),
            genome=str(tmp_path / "genome.fa"),
            gff3=str(tmp_path / "anno.gff3"),
        ),
    )


def _make_stage(tmp_path: Path) -> PropertiesStage:
    config = _make_config(tmp_path)
    fm = FileManager(config)
    fm.ensure_dirs()
    return PropertiesStage(config, fm)


# Two small, real protein sequences for testing
_FASTA_2SEQ = textwrap.dedent("""\
    >gene1
    MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTRDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEECLFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILFLPLPV
    >gene2
    MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRDLYDDDDKDPSSEFQKMAEIAGKDLKDLKDLKDLKDLKDLK
""")

# Sequence with stop codon
_FASTA_STOP = textwrap.dedent("""\
    >gene_stop
    MACDEFGHIKLMNPQRSTVWY*
""")


# ---------------------------------------------------------------------------
# 1. calc_properties with 2 sequences
# ---------------------------------------------------------------------------

def test_calc_properties_two_sequences(tmp_path):
    stage = _make_stage(tmp_path)
    fasta = tmp_path / "test.fa"
    fasta.write_text(_FASTA_2SEQ)

    results = stage.calc_properties(fasta)

    assert len(results) == 2
    assert results[0]["gene_id"] == "gene1"
    assert results[1]["gene_id"] == "gene2"


# ---------------------------------------------------------------------------
# 2. calc_properties returns correct keys
# ---------------------------------------------------------------------------

def test_calc_properties_correct_keys(tmp_path):
    stage = _make_stage(tmp_path)
    fasta = tmp_path / "test.fa"
    fasta.write_text(_FASTA_2SEQ)

    results = stage.calc_properties(fasta)
    row = results[0]

    # scalar columns
    for col in _SCALAR_COLS:
        assert col in row, f"Missing column: {col}"
    # AA composition columns
    for aa in _AA_ORDER:
        assert aa in row, f"Missing amino acid column: {aa}"


# ---------------------------------------------------------------------------
# 3. run raises StageInputError when pep.fa missing
# ---------------------------------------------------------------------------

def test_run_raises_if_pep_missing(tmp_path):
    stage = _make_stage(tmp_path)
    with pytest.raises(StageInputError):
        stage.run()


# ---------------------------------------------------------------------------
# 4. run skips if output exists and force=False
# ---------------------------------------------------------------------------

def test_run_skips_if_output_exists(tmp_path):
    stage = _make_stage(tmp_path)

    stage._pep_fa.write_text(_FASTA_2SEQ)
    stage._output_tsv.write_text("existing content\n")

    stage.run(force=False)

    assert stage._output_tsv.read_text() == "existing content\n"


# ---------------------------------------------------------------------------
# 5. run creates properties.tsv with correct header
# ---------------------------------------------------------------------------

def test_run_creates_tsv_with_correct_header(tmp_path):
    stage = _make_stage(tmp_path)
    stage._pep_fa.write_text(_FASTA_2SEQ)

    stage.run()

    assert stage._output_tsv.exists()
    lines = stage._output_tsv.read_text().splitlines()
    header = lines[0].split("\t")

    for col in _SCALAR_COLS:
        assert col in header, f"Missing scalar column in header: {col}"
    for aa in _AA_ORDER:
        assert aa in header, f"Missing AA column in header: {aa}"


# ---------------------------------------------------------------------------
# 6. calc_properties handles sequences with stop codon (*) correctly
# ---------------------------------------------------------------------------

def test_calc_properties_handles_stop_codon(tmp_path):
    stage = _make_stage(tmp_path)
    fasta = tmp_path / "stop.fa"
    fasta.write_text(_FASTA_STOP)

    results = stage.calc_properties(fasta)

    assert len(results) == 1
    assert results[0]["gene_id"] == "gene_stop"
    assert results[0]["molecular_weight"] > 0
    assert results[0]["isoelectric_point"] > 0


# ---------------------------------------------------------------------------
# 7. run prints online tool hints
# ---------------------------------------------------------------------------

def test_run_prints_online_tool_hints(tmp_path, capsys):
    stage = _make_stage(tmp_path)
    stage._pep_fa.write_text(_FASTA_2SEQ)

    stage.run()

    captured = capsys.readouterr()
    assert "DeepTMHMM" in captured.out
    assert "SignalP" in captured.out
    assert "WoLF PSORT" in captured.out
    assert str(stage._pep_fa) in captured.out


# ---------------------------------------------------------------------------
# 8. calc_properties values are reasonable
# ---------------------------------------------------------------------------

def test_calc_properties_values_reasonable(tmp_path):
    stage = _make_stage(tmp_path)
    fasta = tmp_path / "test.fa"
    fasta.write_text(_FASTA_2SEQ)

    results = stage.calc_properties(fasta)
    row = results[0]

    # Length: positive integer
    assert isinstance(row["length"], int)
    assert row["length"] > 0

    # Molecular weight: positive, reasonable range (Da)
    assert row["molecular_weight"] > 1000

    # Isoelectric point: 0–14
    assert 0 < row["isoelectric_point"] < 14

    # Instability index: typically 0–100 (can exceed 100 for very unstable)
    assert isinstance(row["instability_index"], float)

    # Aliphatic index: positive
    assert row["aliphatic_index"] >= 0

    # GRAVY: typically -2 to +2 for globular proteins
    assert -5 < row["gravy"] < 5

    # Aromaticity: fraction 0–1
    assert 0 <= row["aromaticity"] <= 1

    # AA percentages should sum to ~100.0
    aa_sum = sum(row[aa] for aa in _AA_ORDER)
    assert abs(aa_sum - 100.0) < 0.1


# ---------------------------------------------------------------------------
# 9. run creates correct number of data rows
# ---------------------------------------------------------------------------

def test_run_tsv_row_count(tmp_path):
    stage = _make_stage(tmp_path)
    stage._pep_fa.write_text(_FASTA_2SEQ)

    stage.run()

    lines = [l for l in stage._output_tsv.read_text().splitlines() if l]
    # 1 header + 2 data rows
    assert len(lines) == 3


# ---------------------------------------------------------------------------
# 10. aliphatic index formula correctness
# ---------------------------------------------------------------------------

def test_aliphatic_index_formula(tmp_path):
    """Verify aliphatic index against manual calculation."""
    from gfpipeline.stages.properties import _aliphatic_index

    # All-alanine sequence: AI = 100 * (1.0 + 0 + 0) = 100.0
    aa_pct_ala = {aa: 0.0 for aa in _AA_ORDER}
    aa_pct_ala["A"] = 100.0
    assert abs(_aliphatic_index(aa_pct_ala) - 100.0) < 1e-6

    # All-valine: AI = 100 * (0 + 2.9*1.0 + 0) = 290.0
    aa_pct_val = {aa: 0.0 for aa in _AA_ORDER}
    aa_pct_val["V"] = 100.0
    assert abs(_aliphatic_index(aa_pct_val) - 290.0) < 1e-6
