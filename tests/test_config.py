"""Tests for config/schema.py and config/loader.py.

Covers:
- Valid config parsing
- Missing required fields raises ConfigError with field name
- Default values are filled correctly
- Round-trip serialization (serialize to dict, re-parse, get equivalent object)
"""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest
import yaml

from gfpipeline.config.schema import (
    ToolsConfig,
    DatabasesConfig,
    IdentifyConfig,
    TreeConfig,
    DomainConfig,
    MotifConfig,
    CollinearityConfig,
    GffQcConfig,
    GenomeDbConfig,
    TransConfig,
    PipelineConfig,
)
from gfpipeline.config.loader import load_config
from gfpipeline.core.exceptions import ConfigError


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

MINIMAL_DATABASES = {
    "ref_fa":  "data/ref.fa",
    "genome":  "data/genome.fa",
    "gff3":    "data/anno.gff3",
}

MINIMAL_CONFIG = {
    "project_name": "TestProj",
    "data_dir":     "data/TestProj_data/",
    "result_dir":   "results/TestProj_results/",
    "databases":    MINIMAL_DATABASES,
}


def write_yaml(tmp_path: Path, data: dict, filename: str = "config.yaml") -> Path:
    p = tmp_path / filename
    p.write_text(yaml.dump(data), encoding="utf-8")
    return p


# ---------------------------------------------------------------------------
# Schema unit tests
# ---------------------------------------------------------------------------

class TestToolsConfigDefaults:
    def test_default_muscle(self):
        assert ToolsConfig().muscle == "muscle"

    def test_default_iqtree(self):
        assert ToolsConfig().iqtree == "iqtree2"

    def test_default_mcscanx(self):
        assert ToolsConfig().mcscanx == "MCScanX"

    def test_all_defaults_are_strings(self):
        cfg = ToolsConfig()
        for field in ToolsConfig.model_fields:
            assert isinstance(getattr(cfg, field), str)


class TestDatabasesConfigRequired:
    def test_all_fields_required(self):
        from pydantic import ValidationError
        with pytest.raises(ValidationError):
            DatabasesConfig()  # no fields provided

    def test_valid_databases(self):
        db = DatabasesConfig(**MINIMAL_DATABASES)
        assert db.ref_fa == "data/ref.fa"
        assert db.gff3 == "data/anno.gff3"


class TestIdentifyConfigDefaults:
    def test_hmm_evalue_default(self):
        assert IdentifyConfig().hmm_evalue == 100

    def test_blast_evalue_default(self):
        assert IdentifyConfig().blast_evalue == 1e-5

    def test_blast_threads_default(self):
        assert IdentifyConfig().blast_threads == 10


class TestDomainConfigDefaults:
    def test_genome_cdd_default_none(self):
        assert DomainConfig().genome_cdd is None

    def test_target_domains_default_empty(self):
        assert DomainConfig().target_domains == []


class TestMotifConfigDefaults:
    def test_filter_mode_default(self):
        assert MotifConfig().filter_mode == "any"

    def test_num_motifs_default(self):
        assert MotifConfig().num_motifs == 10


class TestGffQcConfigDefaults:
    def test_rna_bam_default_empty(self):
        assert GffQcConfig().rna_bam == []

    def test_transcript_fa_default_none(self):
        assert GffQcConfig().transcript_fa is None

    def test_output_dir_default(self):
        assert GffQcConfig().output_dir == "data/gff_qc/"


class TestGenomeDbConfigDefaults:
    def test_build_prot_default_true(self):
        assert GenomeDbConfig().build_prot is True

    def test_genome_name_default_none(self):
        assert GenomeDbConfig().genome_name is None

    def test_rep_selection_default(self):
        assert GenomeDbConfig().rep_selection == "longest_cds"


class TestTransConfigDefaults:
    def test_expression_matrix_default_none(self):
        assert TransConfig().expression_matrix is None

    def test_venn_groups_default_empty(self):
        assert TransConfig().venn_groups == []


class TestPipelineConfigDefaults:
    def test_minimal_config_fills_defaults(self):
        cfg = PipelineConfig(**MINIMAL_CONFIG)
        assert cfg.tools == ToolsConfig()
        assert cfg.identify == IdentifyConfig()
        assert cfg.tree == TreeConfig()
        assert cfg.domain == DomainConfig()
        assert cfg.motif == MotifConfig()
        assert cfg.collinearity == CollinearityConfig()
        assert cfg.gff_qc == GffQcConfig()
        assert cfg.genome_db == GenomeDbConfig()
        assert cfg.trans == TransConfig()

    def test_project_name_stored(self):
        cfg = PipelineConfig(**MINIMAL_CONFIG)
        assert cfg.project_name == "TestProj"

    def test_databases_stored(self):
        cfg = PipelineConfig(**MINIMAL_CONFIG)
        assert cfg.databases.ref_fa == "data/ref.fa"


# ---------------------------------------------------------------------------
# load_config tests
# ---------------------------------------------------------------------------

class TestLoadConfigValid:
    def test_load_minimal_config(self, tmp_path):
        p = write_yaml(tmp_path, MINIMAL_CONFIG)
        cfg = load_config(p)
        assert cfg.project_name == "TestProj"
        assert cfg.databases.ref_fa == "data/ref.fa"

    def test_load_config_with_overrides(self, tmp_path):
        data = {**MINIMAL_CONFIG, "identify": {"hmm_evalue": 50, "blast_threads": 4}}
        p = write_yaml(tmp_path, data)
        cfg = load_config(p)
        assert cfg.identify.hmm_evalue == 50
        assert cfg.identify.blast_threads == 4
        # non-overridden default preserved
        assert cfg.identify.blast_evalue == 1e-5

    def test_load_config_accepts_string_path(self, tmp_path):
        p = write_yaml(tmp_path, MINIMAL_CONFIG)
        cfg = load_config(str(p))
        assert isinstance(cfg, PipelineConfig)

    def test_load_config_accepts_path_object(self, tmp_path):
        p = write_yaml(tmp_path, MINIMAL_CONFIG)
        cfg = load_config(p)
        assert isinstance(cfg, PipelineConfig)


class TestLoadConfigMissingFile:
    def test_missing_file_raises_config_error(self, tmp_path):
        with pytest.raises(ConfigError) as exc_info:
            load_config(tmp_path / "nonexistent.yaml")
        assert "not found" in str(exc_info.value).lower()


class TestLoadConfigMissingRequiredFields:
    def test_missing_project_name_raises_config_error(self, tmp_path):
        data = {k: v for k, v in MINIMAL_CONFIG.items() if k != "project_name"}
        p = write_yaml(tmp_path, data)
        with pytest.raises(ConfigError) as exc_info:
            load_config(p)
        assert "project_name" in str(exc_info.value)

    def test_missing_data_dir_raises_config_error(self, tmp_path):
        data = {k: v for k, v in MINIMAL_CONFIG.items() if k != "data_dir"}
        p = write_yaml(tmp_path, data)
        with pytest.raises(ConfigError) as exc_info:
            load_config(p)
        assert "data_dir" in str(exc_info.value)

    def test_missing_result_dir_raises_config_error(self, tmp_path):
        data = {k: v for k, v in MINIMAL_CONFIG.items() if k != "result_dir"}
        p = write_yaml(tmp_path, data)
        with pytest.raises(ConfigError) as exc_info:
            load_config(p)
        assert "result_dir" in str(exc_info.value)

    def test_missing_databases_raises_config_error(self, tmp_path):
        data = {k: v for k, v in MINIMAL_CONFIG.items() if k != "databases"}
        p = write_yaml(tmp_path, data)
        with pytest.raises(ConfigError) as exc_info:
            load_config(p)
        assert "databases" in str(exc_info.value)

    def test_missing_databases_ref_fa_raises_config_error(self, tmp_path):
        db = {k: v for k, v in MINIMAL_DATABASES.items() if k != "ref_fa"}
        data = {**MINIMAL_CONFIG, "databases": db}
        p = write_yaml(tmp_path, data)
        with pytest.raises(ConfigError) as exc_info:
            load_config(p)
        assert "ref_fa" in str(exc_info.value)

    def test_missing_databases_gff3_raises_config_error(self, tmp_path):
        db = {k: v for k, v in MINIMAL_DATABASES.items() if k != "gff3"}
        data = {**MINIMAL_CONFIG, "databases": db}
        p = write_yaml(tmp_path, data)
        with pytest.raises(ConfigError) as exc_info:
            load_config(p)
        assert "gff3" in str(exc_info.value)

    def test_error_is_config_error_subclass(self, tmp_path):
        data = {k: v for k, v in MINIMAL_CONFIG.items() if k != "project_name"}
        p = write_yaml(tmp_path, data)
        with pytest.raises(ConfigError):
            load_config(p)

    def test_invalid_yaml_raises_config_error(self, tmp_path):
        p = tmp_path / "bad.yaml"
        p.write_text("key: [unclosed", encoding="utf-8")
        with pytest.raises(ConfigError):
            load_config(p)


# ---------------------------------------------------------------------------
# Round-trip serialization tests
# ---------------------------------------------------------------------------

class TestRoundTrip:
    def test_roundtrip_minimal_config(self):
        cfg = PipelineConfig(**MINIMAL_CONFIG)
        data = cfg.model_dump()
        cfg2 = PipelineConfig(**data)
        assert cfg == cfg2

    def test_roundtrip_with_custom_tools(self):
        data = {**MINIMAL_CONFIG, "tools": {"muscle": "/usr/local/bin/muscle"}}
        cfg = PipelineConfig(**data)
        cfg2 = PipelineConfig(**cfg.model_dump())
        assert cfg2.tools.muscle == "/usr/local/bin/muscle"
        assert cfg == cfg2

    def test_roundtrip_with_optional_fields(self):
        data = {
            **MINIMAL_CONFIG,
            "domain": {"genome_cdd": "data/genome.cdd.txt", "target_domains": ["cd12345"]},
            "trans": {"expression_matrix": "data/expr.tsv", "logfc_threshold": 1.5},
        }
        cfg = PipelineConfig(**data)
        cfg2 = PipelineConfig(**cfg.model_dump())
        assert cfg == cfg2

    def test_roundtrip_via_yaml_file(self, tmp_path):
        cfg = PipelineConfig(**MINIMAL_CONFIG)
        # serialize to YAML and reload
        raw = cfg.model_dump()
        p = write_yaml(tmp_path, raw)
        cfg2 = load_config(p)
        assert cfg == cfg2

    def test_roundtrip_preserves_list_fields(self):
        data = {
            **MINIMAL_CONFIG,
            "domain": {"target_domains": ["cd001", "cd002"]},
            "gff_qc": {"rna_bam": ["a.bam", "b.bam"], "extra_gff": ["extra.gff3"]},
        }
        cfg = PipelineConfig(**data)
        cfg2 = PipelineConfig(**cfg.model_dump())
        assert cfg2.domain.target_domains == ["cd001", "cd002"]
        assert cfg2.gff_qc.rna_bam == ["a.bam", "b.bam"]
        assert cfg == cfg2
