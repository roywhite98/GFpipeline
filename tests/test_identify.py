"""Unit tests for stages/identify.py — IdentifyStage."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from gfpipeline.config.schema import DatabasesConfig, GenomeDbConfig, PipelineConfig
from gfpipeline.core.exceptions import StageInputError
from gfpipeline.core.file_manager import FileManager
from gfpipeline.core.runner import ToolRunner
from gfpipeline.stages.identify import IdentifyStage


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def make_config(tmp_path: Path) -> PipelineConfig:
    data_dir = tmp_path / "data"
    result_dir = tmp_path / "results"
    data_dir.mkdir()
    result_dir.mkdir()
    return PipelineConfig(
        project_name="TEST",
        data_dir=str(data_dir),
        result_dir=str(result_dir),
        databases=DatabasesConfig(
            ref_fa=str(tmp_path / "ref.fa"),
            genome=str(tmp_path / "genome.fa"),
            gff3=str(tmp_path / "anno.gff3"),
        ),
        genome_db=GenomeDbConfig(rep_index_dir=str(tmp_path / "rep_index")),
    )


def make_stage(tmp_path: Path) -> tuple[IdentifyStage, PipelineConfig]:
    config = make_config(tmp_path)
    runner = ToolRunner(dry_run=True)
    fm = FileManager(config)
    stage = IdentifyStage(config, runner, fm)
    return stage, config


# ---------------------------------------------------------------------------
# merge_ids
# ---------------------------------------------------------------------------

class TestMergeIds:
    def test_union_of_both_lists(self, tmp_path):
        stage, _ = make_stage(tmp_path)
        result = stage.merge_ids(["A", "B"], ["C", "D"])
        assert set(result) == {"A", "B", "C", "D"}

    def test_deduplicates_overlapping_ids(self, tmp_path):
        stage, _ = make_stage(tmp_path)
        result = stage.merge_ids(["A", "B"], ["B", "C"])
        assert result.count("B") == 1

    def test_result_is_sorted(self, tmp_path):
        stage, _ = make_stage(tmp_path)
        result = stage.merge_ids(["C", "A"], ["B"])
        assert result == sorted(result)

    def test_empty_hmm_ids(self, tmp_path):
        stage, _ = make_stage(tmp_path)
        result = stage.merge_ids([], ["X", "Y"])
        assert result == ["X", "Y"]

    def test_empty_blast_ids(self, tmp_path):
        stage, _ = make_stage(tmp_path)
        result = stage.merge_ids(["X", "Y"], [])
        assert result == ["X", "Y"]

    def test_both_empty(self, tmp_path):
        stage, _ = make_stage(tmp_path)
        result = stage.merge_ids([], [])
        assert result == []

    def test_result_length_le_sum_of_inputs(self, tmp_path):
        stage, _ = make_stage(tmp_path)
        hmm = ["A", "B", "C"]
        blast = ["B", "C", "D"]
        result = stage.merge_ids(hmm, blast)
        assert len(result) <= len(hmm) + len(blast)

    def test_all_ids_present(self, tmp_path):
        stage, _ = make_stage(tmp_path)
        hmm = ["Os01t0001-01", "Os02t0002-01"]
        blast = ["Os03t0003-01", "Os01t0001-01"]
        result = stage.merge_ids(hmm, blast)
        for id_ in hmm + blast:
            assert id_ in result


# ---------------------------------------------------------------------------
# Input validation (StageInputError)
# ---------------------------------------------------------------------------

class TestRunInputValidation:
    def test_raises_when_neither_ref_nor_hmm_exists(self, tmp_path):
        stage, _ = make_stage(tmp_path)
        with pytest.raises(StageInputError):
            stage.run()

    def test_no_error_when_hmm_exists(self, tmp_path):
        stage, config = make_stage(tmp_path)
        # hmm is in same dir as ref_fa (data_dir)
        hmm_path = stage._hmm
        hmm_path.parent.mkdir(parents=True, exist_ok=True)
        hmm_path.write_text("dummy hmm")
        # Also need rep_pep to exist
        stage.fm.rep_pep.parent.mkdir(parents=True, exist_ok=True)
        stage.fm.rep_pep.write_text(">seq1\nMAAA\n")
        with (
            patch.object(stage, "hmm_search", return_value=[]),
            patch.object(stage, "blast_search", return_value=[]),
            patch.object(stage, "extract_sequences"),
        ):
            # Should not raise StageInputError
            stage.run()

    def test_no_error_when_ref_fa_exists(self, tmp_path):
        stage, config = make_stage(tmp_path)
        ref_fa = Path(config.databases.ref_fa)
        ref_fa.parent.mkdir(parents=True, exist_ok=True)
        ref_fa.write_text(">seq1\nMAAA\n")
        # Also need rep_pep to exist
        stage.fm.rep_pep.parent.mkdir(parents=True, exist_ok=True)
        stage.fm.rep_pep.write_text(">seq1\nMAAA\n")
        with (
            patch.object(stage, "build_hmm"),
            patch.object(stage, "hmm_search", return_value=[]),
            patch.object(stage, "blast_search", return_value=[]),
            patch.object(stage, "extract_sequences"),
        ):
            stage.run()


# ---------------------------------------------------------------------------
# Path helpers
# ---------------------------------------------------------------------------

class TestPathHelpers:
    def test_ref_fa_path(self, tmp_path):
        stage, config = make_stage(tmp_path)
        assert stage._ref_fa == Path(config.databases.ref_fa)

    def test_hmm_path(self, tmp_path):
        stage, config = make_stage(tmp_path)
        # ref.fa → stem "ref" → ref.hmm (stem derived from ref_fa filename)
        assert stage._hmm == Path(config.databases.ref_fa).parent / "ref.hmm"

    def test_hmm_path_ref_stem_convention(self, tmp_path):
        # wrky.ref.fa → stem "wrky" → wrky.hmm
        config = make_config(tmp_path)
        config = config.model_copy(update={
            "databases": config.databases.model_copy(update={
                "ref_fa": str(tmp_path / "wrky.ref.fa")
            })
        })
        runner = ToolRunner(dry_run=True)
        stage = IdentifyStage(config, runner, FileManager(config))
        assert stage._hmm == tmp_path / "wrky.hmm"
        assert stage._hmmemit_out == tmp_path / "wrky.hmmemit.out"

    def test_candidates_pep_fa_path(self, tmp_path):
        stage, config = make_stage(tmp_path)
        expected = Path(config.result_dir) / "TEST.identify.candidates.pep.fa"
        assert stage._candidates_pep_fa == expected

    def test_candidates_gene_idlist_path(self, tmp_path):
        stage, config = make_stage(tmp_path)
        expected = Path(config.result_dir) / "TEST.identify.candidates.gene.idlist"
        assert stage._candidates_gene_idlist == expected


# ---------------------------------------------------------------------------
# build_hmm (dry_run — commands not executed, files not created)
# ---------------------------------------------------------------------------

class TestBuildHmm:
    def test_build_hmm_calls_muscle_and_hmmbuild(self, tmp_path):
        stage, config = make_stage(tmp_path)
        ref_fa = Path(config.data_dir) / "TEST.ref.fa"
        ref_fa.write_text(">seq1\nMAAA\n")

        called_cmds = []
        original_run = stage.runner.run

        def capture_run(cmd, **kwargs):
            called_cmds.append(cmd[0])
            return original_run(cmd, **kwargs)

        stage.runner.run = capture_run
        stage.build_hmm()

        assert "muscle" in called_cmds
        assert "hmmbuild" in called_cmds

    def test_build_hmm_skips_muscle_if_afa_exists(self, tmp_path):
        stage, config = make_stage(tmp_path)
        ref_fa = Path(config.data_dir) / "TEST.ref.fa"
        ref_fa.write_text(">seq1\nMAAA\n")
        # Pre-create the afa file
        afa = Path(config.result_dir) / "TEST.identify.ref.afa"
        afa.write_text("dummy afa")

        called_cmds = []
        original_run = stage.runner.run

        def capture_run(cmd, **kwargs):
            called_cmds.append(cmd[0])
            return original_run(cmd, **kwargs)

        stage.runner.run = capture_run
        stage.build_hmm()

        assert "muscle" not in called_cmds


# ---------------------------------------------------------------------------
# extract_sequences
# ---------------------------------------------------------------------------

FASTA_PEP = """\
>Os01t0001-01 protein1
MAAAKLL
>Os02t0002-01 protein2
MTTTTLL
"""

FASTA_CDS = """\
>Os01t0001-01 cds1
ATGAAAGCG
>Os02t0002-01 cds2
ATGACAACA
"""

FASTA_GENOME = """\
>Os01g0001 gene1
AAATTTGGG
>Os02g0002 gene2
CCCGGGAAA
"""


class TestExtractSequences:
    def test_extract_sequences_creates_output_files(self, tmp_path):
        stage, config = make_stage(tmp_path)

        # Write real FASTA databases to rep_pep/rep_cds locations
        stage.fm.rep_pep.parent.mkdir(parents=True, exist_ok=True)
        stage.fm.rep_pep.write_text(FASTA_PEP)
        stage.fm.rep_cds.write_text(FASTA_CDS)
        Path(config.databases.genome).write_text(FASTA_GENOME)

        ids = ["Os01t0001-01", "Os02t0002-01"]
        gene_ids = ["Os01g0001", "Os02g0002"]

        stage.extract_sequences(ids, gene_ids)

        assert stage._candidates_pep_fa.exists()
        assert stage._candidates_cds_fa.exists()
        assert stage._candidates_genome_fa.exists()

    def test_extract_sequences_correct_content(self, tmp_path):
        stage, config = make_stage(tmp_path)

        stage.fm.rep_pep.parent.mkdir(parents=True, exist_ok=True)
        stage.fm.rep_pep.write_text(FASTA_PEP)
        stage.fm.rep_cds.write_text(FASTA_CDS)
        Path(config.databases.genome).write_text(FASTA_GENOME)

        stage.extract_sequences(["Os01t0001-01"], ["Os01g0001"])

        pep_content = stage._candidates_pep_fa.read_text()
        assert "Os01t0001-01" in pep_content
        assert "Os02t0002-01" not in pep_content
