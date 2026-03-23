"""Unit tests for core/tool_checker.py."""

import sys
from unittest.mock import patch

import pytest

from gfpipeline.config.schema import DatabasesConfig, PipelineConfig
from gfpipeline.core.exceptions import ToolNotFoundError
from gfpipeline.core.tool_checker import STAGE_TOOLS, check_tools, is_executable


def make_config(**kwargs) -> PipelineConfig:
    defaults = dict(
        project_name="TEST",
        data_dir="data/",
        result_dir="results/",
        databases=DatabasesConfig(
            ref_fa="ref.fa",
            genome="genome.fa",
            gff3="anno.gff3",
        ),
    )
    defaults.update(kwargs)
    return PipelineConfig(**defaults)


# ---------------------------------------------------------------------------
# is_executable
# ---------------------------------------------------------------------------

class TestIsExecutable:
    def test_python_is_executable(self):
        # python interpreter is always available
        assert is_executable(sys.executable) is True

    def test_nonexistent_tool_returns_false(self):
        assert is_executable("__nonexistent_tool_xyz__") is False

    def test_path_to_nonexistent_file_returns_false(self, tmp_path):
        assert is_executable(str(tmp_path / "no_such_tool")) is False

    def test_path_to_non_executable_file_returns_false(self, tmp_path):
        f = tmp_path / "script.sh"
        f.write_text("#!/bin/sh\necho hi")
        # not chmod +x
        assert is_executable(str(f)) is False

    def test_path_to_executable_file_returns_true(self, tmp_path):
        f = tmp_path / "mytool"
        f.write_text("#!/bin/sh\necho hi")
        f.chmod(0o755)
        assert is_executable(str(f)) is True


# ---------------------------------------------------------------------------
# STAGE_TOOLS mapping
# ---------------------------------------------------------------------------

class TestStageTools:
    def test_identify_has_required_tools(self):
        tools = STAGE_TOOLS["identify"]
        for t in ["muscle", "hmmbuild", "hmmsearch", "hmmemit", "blastp"]:
            assert t in tools

    def test_tree_has_required_tools(self):
        tools = STAGE_TOOLS["tree"]
        for t in ["muscle", "trimal", "iqtree2"]:
            assert t in tools

    def test_domain_has_no_tools(self):
        assert STAGE_TOOLS["domain"] == []

    def test_collinearity_has_blastp(self):
        assert "blastp" in STAGE_TOOLS["collinearity"]

    def test_genome_db_blast_has_makeblastdb(self):
        assert "makeblastdb" in STAGE_TOOLS["genome-db blast"]

    def test_genome_db_gff_qc_has_samtools(self):
        assert "samtools" in STAGE_TOOLS["genome-db gff-qc"]

    def test_genome_db_rep_index_is_empty(self):
        assert STAGE_TOOLS["genome-db rep-index"] == []


# ---------------------------------------------------------------------------
# check_tools
# ---------------------------------------------------------------------------

class TestCheckTools:
    def test_domain_stage_passes_with_no_tools(self):
        config = make_config()
        # domain stage has no tools, should never raise
        check_tools("domain", config)

    def test_genome_db_rep_index_passes_with_no_tools(self):
        config = make_config()
        check_tools("genome-db rep-index", config)

    def test_raises_tool_not_found_for_missing_tool(self):
        config = make_config()
        # Override muscle to a nonexistent path
        config.tools.muscle = "__nonexistent_muscle__"
        with pytest.raises(ToolNotFoundError) as exc_info:
            check_tools("identify", config)
        assert "__nonexistent_muscle__" in str(exc_info.value)

    def test_error_message_contains_stage_name(self):
        config = make_config()
        config.tools.muscle = "__nonexistent_muscle__"
        with pytest.raises(ToolNotFoundError) as exc_info:
            check_tools("identify", config)
        assert "identify" in str(exc_info.value)

    def test_collinearity_mcscanx_adds_mcscanx_check(self):
        from gfpipeline.config.schema import CollinearityConfig
        config = make_config()
        config.collinearity = CollinearityConfig(tool="mcscanx")
        config.tools.mcscanx = "__nonexistent_mcscanx__"
        # blastp is also checked; mock it as available
        with patch("gfpipeline.core.tool_checker.is_executable") as mock_exec:
            mock_exec.side_effect = lambda t: t != "__nonexistent_mcscanx__"
            with pytest.raises(ToolNotFoundError) as exc_info:
                check_tools("collinearity", config)
            assert "__nonexistent_mcscanx__" in str(exc_info.value)

    def test_collinearity_jcvi_adds_jcvi_check(self):
        from gfpipeline.config.schema import CollinearityConfig
        config = make_config()
        config.collinearity = CollinearityConfig(tool="jcvi")
        with patch("gfpipeline.core.tool_checker.is_executable") as mock_exec:
            # All tools pass except jcvi
            mock_exec.side_effect = lambda t: t != "jcvi"
            with pytest.raises(ToolNotFoundError) as exc_info:
                check_tools("collinearity", config)
            assert "jcvi" in str(exc_info.value)

    def test_gff_qc_adds_minimap2_when_transcript_fa_configured(self):
        from gfpipeline.config.schema import GffQcConfig
        config = make_config()
        config.gff_qc = GffQcConfig(transcript_fa="transcripts.fa")
        config.tools.minimap2 = "__nonexistent_minimap2__"
        with patch("gfpipeline.core.tool_checker.is_executable") as mock_exec:
            mock_exec.side_effect = lambda t: t != "__nonexistent_minimap2__"
            with pytest.raises(ToolNotFoundError) as exc_info:
                check_tools("genome-db gff-qc", config)
            assert "__nonexistent_minimap2__" in str(exc_info.value)

    def test_gff_qc_adds_stringtie_when_rna_bam_configured(self):
        from gfpipeline.config.schema import GffQcConfig
        config = make_config()
        config.gff_qc = GffQcConfig(rna_bam=["sample.bam"])
        config.tools.stringtie = "__nonexistent_stringtie__"
        with patch("gfpipeline.core.tool_checker.is_executable") as mock_exec:
            mock_exec.side_effect = lambda t: t != "__nonexistent_stringtie__"
            with pytest.raises(ToolNotFoundError) as exc_info:
                check_tools("genome-db gff-qc", config)
            assert "__nonexistent_stringtie__" in str(exc_info.value)

    def test_gff_qc_no_optional_tools_when_not_configured(self):
        from gfpipeline.config.schema import GffQcConfig
        config = make_config()
        config.gff_qc = GffQcConfig()  # no transcript_fa, no rna_bam
        # Only samtools is required; mock it as available
        with patch("gfpipeline.core.tool_checker.is_executable", return_value=True):
            check_tools("genome-db gff-qc", config)  # should not raise
