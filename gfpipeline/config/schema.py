"""Pydantic configuration schema for gfpipeline."""

from __future__ import annotations

from typing import Optional

from pydantic import BaseModel


class ToolsConfig(BaseModel):
    muscle:      str = "muscle"
    hmmbuild:    str = "hmmbuild"
    hmmsearch:   str = "hmmsearch"
    hmmemit:     str = "hmmemit"
    blastp:      str = "blastp"
    trimal:      str = "trimal"
    iqtree:      str = "iqtree2"
    meme:        str = "meme"
    fimo:        str = "fimo"
    mcscanx:     str = "MCScanX"
    makeblastdb: str = "makeblastdb"
    samtools:    str = "samtools"
    minimap2:    str = "minimap2"
    stringtie:   str = "stringtie"


class DatabasesConfig(BaseModel):
    ref_fa:  str  # 目标基因家族参考蛋白序列 FASTA
    genome:  str  # 全基因组序列 FASTA
    gff3:    str  # GFF3 注释文件


class IdentifyConfig(BaseModel):
    hmm_evalue:    float = 100
    blast_evalue:  float = 1e-5
    blast_threads: int   = 10


class TreeConfig(BaseModel):
    iqtree_bootstrap: int = 1000


class DomainConfig(BaseModel):
    evalue:         float            = 0.01
    maxhit:         int              = 250
    genome_cdd:     Optional[str]    = None
    target_domains: list[str]        = []


class MotifConfig(BaseModel):
    num_motifs:      int = 10
    min_width:       int = 6
    max_width:       int = 50
    filter_mode:     str = "any"   # any | all | min_count
    min_motif_count: int = 1


class CollinearityConfig(BaseModel):
    tool:          str = "jcvi"   # jcvi | mcscanx
    blast_threads: int = 10


class GffQcConfig(BaseModel):
    output_dir:    str           = "data/gff_qc/"
    min_cds_len:   int           = 150
    edge_distance: int           = 500
    rna_bam:       list[str]     = []
    transcript_fa: Optional[str] = None
    extra_gff:     list[str]     = []


class GenomeDbConfig(BaseModel):
    build_prot:    bool          = True
    build_nucl:    bool          = True
    index_dir:     str           = "data/genome_index/"
    rep_index_dir: str           = "data/rep_index/"
    genome_name:   Optional[str] = None   # 默认使用 project_name
    rep_selection: str           = "longest_cds"  # longest_cds | longest_mrna


class TransConfig(BaseModel):
    expression_matrix: Optional[str]   = None
    logfc_threshold:   Optional[float] = None
    pvalue_threshold:  Optional[float] = None
    venn_groups:       list[str]       = []


class PipelineConfig(BaseModel):
    project_name: str
    data_dir:     str
    result_dir:   str
    databases:    DatabasesConfig
    tools:        ToolsConfig        = ToolsConfig()
    identify:     IdentifyConfig     = IdentifyConfig()
    tree:         TreeConfig         = TreeConfig()
    domain:       DomainConfig       = DomainConfig()
    motif:        MotifConfig        = MotifConfig()
    collinearity: CollinearityConfig = CollinearityConfig()
    gff_qc:       GffQcConfig        = GffQcConfig()
    genome_db:    GenomeDbConfig     = GenomeDbConfig()
    trans:        TransConfig        = TransConfig()
