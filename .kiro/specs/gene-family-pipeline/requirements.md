# 需求文档

## 简介

本功能将现有的基因家族分析 Jupyter Notebook 工作流（IDENTIFY.ipynb、Tree.ipynb、WRF.ipynb、PROPERTIES.ipynb、TRANS.ipynb）重构为一个结构化的 Python CLI 应用（pipeline）。该 CLI 应用支持按阶段独立运行或一键执行全流程，通过配置文件驱动项目参数，结果统一输出到 `results/{Proj_Name}_results/` 目录，从而提升分析流程的可重复性、可维护性和自动化程度。

流程涵盖以下阶段：基因组数据库准备（genome-db，前置阶段）、基因家族鉴定（identify）、系统发育树构建（tree）、保守结构域分析（domain）、结构域筛选（domain-filter）、Motif 发现与筛选（motif）、基因组共线性分析（collinearity）、基因属性分析（properties）、转录组分析（trans）。

genome-db 阶段包含四个子步骤：BLAST 数据库构建（blast）、GFF 文件质量检查与修正（gff-qc）、基因序列索引构建（gene-index）、代表基因索引构建（rep-index）。这四个子步骤可作为 identify 阶段的前置准备统一运行，也可独立调用以处理新的基因组数据。

---

## 词汇表

- **Pipeline**：本 CLI 应用的整体分析流程，由多个阶段（Stage）组成。
- **Stage**：流程中的一个独立分析阶段，可单独运行，也可作为全流程的一部分运行。
- **Proj_Name**：项目名称，用于命名数据目录、结果目录及中间文件。
- **Config**：配置文件（YAML 格式），存储项目名称、数据库路径、工具路径及各阶段参数。
- **HMM_Builder**：负责执行多序列比对（muscle）和 HMM profile 构建（hmmbuild）的模块。
- **HMM_Searcher**：负责执行 hmmsearch 并提取候选成员 ID 的模块。
- **BLAST_Searcher**：负责执行 hmmemit + blastp 并提取候选成员 ID 的模块。
- **ID_Merger**：负责合并 HMM 和 BLAST 结果、去重、转换 transcript ID 为 gene ID 的模块。
- **Sequence_Extractor**：负责根据 ID 列表从 FASTA 数据库提取 CDS、蛋白、基因组序列的模块。
- **Tree_Builder**：负责执行多序列比对、trimAl 修剪和 IQ-TREE 建树的模块。
- **Domain_Analyzer**：负责提交序列到 NCBI Batch CD-Search 并解析结果的模块。
- **Domain_Filter**：负责解析 CDD 结果、提取目标结构域 accession、从全基因组 CDD 结果中筛选候选基因的模块。
- **Motif_Analyzer**：负责调用 MEME 发现 motif、调用 FIMO 扫描全基因组蛋白、解析结果并筛选候选基因的模块。
- **Collinearity_Analyzer**：负责执行全基因组自比对 blastp、调用 MCScanX 或 JCVI 进行共线性分析、提取包含目标基因的共线性块的模块。
- **Summary_Reporter**：负责在全流程结束后汇总各阶段产出文件路径和统计数字的模块。
- **Result_Dir**：结果输出目录，路径为 `results/{Proj_Name}_results/`。
- **Data_Dir**：项目数据目录，路径为 `data/{Proj_Name}_data/`。
- **Transcript_ID**：转录本级别的序列标识符（如 `Os01t0936800-01`）。
- **Gene_ID**：基因级别的标识符（如 `Os01g0936800`），由 Transcript_ID 转换而来。
- **CDD**：NCBI 保守结构域数据库（Conserved Domain Database）。
- **Domain_Accession**：CDD 中结构域的唯一标识符（如 `cd12345`）。
- **Superfamily_Accession**：CDD 中超家族的唯一标识符（如 `cl12345`）。
- **Motif**：由 MEME 发现的蛋白序列保守模式。
- **Collinearity_Block**：基因组共线性分析中识别出的保守基因顺序片段。
- **Genome_DB**：基因组数据库准备阶段（genome-db），包含 BLAST 数据库构建、GFF 质检修正、基因序列索引和代表基因索引四个子步骤，为后续 identify 阶段提供标准化数据库。
- **BLAST_DB_Builder**：负责调用 makeblastdb 为蛋白或核酸 FASTA 文件构建 BLAST 数据库的模块。
- **GFF_QC**：负责对 GFF3 文件进行格式规范性检查、基因模型完整性检查、疑似不完整基因标记及修正的模块。
- **Gene_Index_Builder**：负责为基因组、CDS、蛋白 FASTA 文件建立序列索引，并构建基因 ID 与转录本 ID 映射表的模块。
- **Rep_Index_Builder**：负责从每个基因的所有转录本中选取代表转录本，并提取代表转录本序列的模块。
- **Representative_Transcript**：每个基因中 CDS 最长的转录本（若 CDS 等长则取 transcript ID 字典序最小者），用于代表该基因参与后续分析。
- **Truncated_Gene**：疑似不完整的基因模型，判定条件包括：位于 scaffold/contig 末端 `gff_qc.edge_distance` bp 范围内、CDS 长度短于 `gff_qc.min_cds_len`、或缺少 start/stop codon。
- **GFF3**：基因组特征注释文件格式（Generic Feature Format Version 3），包含基因、转录本、CDS、exon 等特征的坐标和属性信息。
- **Index_Dir**：基因序列索引输出目录，路径由配置参数 `genome_db.index_dir` 指定，默认为 `data/genome_index/`。
- **Rep_Index_Dir**：代表基因索引输出目录，路径由配置参数 `genome_db.rep_index_dir` 指定，默认为 `data/rep_index/`。
- **Genome_Name**：用于命名 genome-db 阶段输出文件的标识符，由配置参数 `genome_db.genome_name` 指定，默认使用 `project_name`。

---

## 需求

### 需求 1：CLI 入口与全局参数

**用户故事：** 作为生物信息学研究人员，我希望通过统一的命令行入口运行分析流程，以便灵活控制运行哪些阶段。

#### 验收标准

1. THE Pipeline SHALL 提供名为 `gfpipeline` 的命令行入口，支持子命令 `run`、`genome-db`、`identify`、`tree`、`domain`、`domain-filter`、`motif`、`collinearity`、`properties`、`trans`。
2. THE Pipeline SHALL 支持 `genome-db` 子命令下的子子命令 `blast`、`gff-qc`、`gene-index`、`rep-index`、`query`。
3. WHEN 用户执行 `gfpipeline genome-db --config <config_file>` 且未指定子子命令时，THE Pipeline SHALL 按顺序依次执行 blast、gff-qc、gene-index、rep-index 四个子步骤。
4. WHEN 用户执行 `gfpipeline run --config <config_file>` 时，THE Pipeline SHALL 按顺序依次执行 identify、tree、domain、domain-filter、motif、collinearity 六个阶段；WHEN `databases.blast_db` 对应的 BLAST 数据库文件不存在时，THE Pipeline SHALL 在执行 identify 之前自动先执行 `genome-db blast` 步骤。
5. WHEN 用户执行 `gfpipeline <stage> --config <config_file>` 时，THE Pipeline SHALL 仅执行指定阶段。
6. IF 用户未提供 `--config` 参数，THEN THE Pipeline SHALL 在当前目录查找名为 `config.yaml` 的默认配置文件；若仍未找到，则输出错误信息并以非零状态码退出。
7. THE Pipeline SHALL 支持 `--dry-run` 标志，WHEN 该标志被设置时，THE Pipeline SHALL 打印将要执行的命令但不实际执行。
8. THE Pipeline SHALL 支持 `--verbose` 标志，WHEN 该标志被设置时，THE Pipeline SHALL 输出每个步骤的详细日志（包括执行的命令、输入输出文件路径）。
9. WHEN 任意阶段执行失败时，THE Pipeline SHALL 输出包含阶段名称和失败原因的错误信息，并以非零状态码退出。

---

### 需求 2：配置文件解析

**用户故事：** 作为研究人员，我希望通过 YAML 配置文件集中管理项目参数，以便在不修改代码的情况下切换项目。

#### 验收标准

1. THE Config_Parser SHALL 解析 YAML 格式的配置文件，并将其映射为结构化的配置对象。
2. THE Config_Parser SHALL 支持以下必填字段：`project_name`、`data_dir`、`result_dir`、`databases.pep`、`databases.cds`、`databases.genome`、`databases.blast_db`、`databases.gff3`（GFF3 注释文件路径，用于共线性分析）。
3. THE Config_Parser SHALL 支持以下工具路径可选字段及其默认值：`tools.muscle`（默认 `muscle`）、`tools.hmmbuild`（默认 `hmmbuild`）、`tools.hmmsearch`（默认 `hmmsearch`）、`tools.hmmemit`（默认 `hmmemit`）、`tools.blastp`（默认 `blastp`）、`tools.trimal`（默认 `trimal`）、`tools.iqtree`（默认 `iqtree2`）、`tools.meme`（默认 `meme`）、`tools.fimo`（默认 `fimo`）、`tools.mcscanx`（默认 `MCScanX`）、`tools.makeblastdb`（默认 `makeblastdb`）、`tools.samtools`（默认 `samtools`）、`tools.minimap2`（默认 `minimap2`，可选，用于 gff-qc 转录本比对）、`tools.stringtie`（默认 `stringtie`，可选，用于 gff-qc 转录本重预测）。
4. THE Config_Parser SHALL 支持各阶段参数覆盖，包括：`identify.hmm_evalue`（默认 `100`）、`identify.blast_evalue`（默认 `1e-5`）、`identify.blast_threads`（默认 `10`）、`tree.iqtree_bootstrap`（默认 `1000`）、`domain.evalue`（默认 `0.01`）、`domain.maxhit`（默认 `250`）、`domain.genome_cdd`（全基因组 CDD 结果路径，可选）、`domain.target_domains`（目标结构域 accession 列表，可选）、`motif.num_motifs`（默认 `10`）、`motif.min_width`（默认 `6`）、`motif.max_width`（默认 `50`）、`motif.filter_mode`（默认 `any`，可选 `any`/`all`/`min_count`）、`motif.min_motif_count`（默认 `1`）、`collinearity.tool`（默认 `jcvi`，可选 `jcvi`/`mcscanx`）、`collinearity.blast_threads`（默认 `10`）、`genome_db.build_prot`（默认 `true`）、`genome_db.build_nucl`（默认 `true`）、`genome_db.index_dir`（默认 `data/genome_index/`）、`genome_db.rep_index_dir`（默认 `data/rep_index/`）、`genome_db.genome_name`（默认使用 `project_name`）、`genome_db.rep_selection`（默认 `longest_cds`，可选 `longest_cds`/`longest_mrna`）、`gff_qc.output_dir`（默认 `data/gff_qc/`）、`gff_qc.min_cds_len`（默认 `150`）、`gff_qc.edge_distance`（默认 `500`）、`gff_qc.rna_bam`（可选，RNA-seq BAM 文件路径列表）、`gff_qc.transcript_fa`（可选，转录本 FASTA 路径）、`gff_qc.extra_gff`（可选，外部 GFF3 路径列表）。
5. IF 配置文件中缺少任意必填字段，THEN THE Config_Parser SHALL 输出包含缺失字段名称的错误信息，并以非零状态码退出。
6. FOR ALL 有效的配置文件，THE Config_Parser 解析后再序列化再解析 SHALL 产生等价的配置对象（往返属性）。

---

### 需求 3：阶段一——基因家族鉴定（identify）

**用户故事：** 作为研究人员，我希望通过 HMM 搜索和 BLAST 搜索联合鉴定候选基因家族成员，以便获得全面的候选集合。

#### 验收标准

1. WHEN identify 阶段启动时，THE Pipeline SHALL 创建 `Data_Dir` 和 `Result_Dir` 目录（若已存在则跳过）。
2. WHEN `Data_Dir/{Proj_Name}.ref.fa` 存在且 `Data_Dir/{Proj_Name}.hmm` 不存在时，THE HMM_Builder SHALL 调用 muscle 对参考蛋白序列进行多序列比对，输出比对文件到 `Result_Dir/{Proj_Name}.identify.ref.afa`。
3. WHEN 多序列比对完成后，THE HMM_Builder SHALL 调用 hmmbuild 构建 HMM profile，输出到 `Data_Dir/{Proj_Name}.hmm`。
4. WHEN `Data_Dir/{Proj_Name}.hmm` 存在时，THE HMM_Searcher SHALL 调用 hmmsearch 在蛋白数据库中搜索，输出原始结果到 `Result_Dir/{Proj_Name}.identify.hmm.out`。
5. WHEN hmmsearch 完成后，THE HMM_Searcher SHALL 从搜索结果中提取 e-value 低于 `identify.hmm_evalue` 阈值的候选 ID，输出到 `Result_Dir/{Proj_Name}.identify.hmm.idlist`。
6. WHEN HMM 文件存在时，THE BLAST_Searcher SHALL 调用 hmmemit 生成 consensus 序列，输出到 `Data_Dir/{Proj_Name}.hmmemit.out`。
7. WHEN consensus 序列生成后，THE BLAST_Searcher SHALL 调用 blastp 在 BLAST 数据库中搜索，e-value 阈值为 `identify.blast_evalue`，输出原始结果到 `Result_Dir/{Proj_Name}.identify.blast.out`。
8. WHEN blastp 完成后，THE BLAST_Searcher SHALL 从 tabular 格式输出中提取第二列（subject ID），去重后输出到 `Result_Dir/{Proj_Name}.identify.blast.idlist`。
9. WHEN HMM 和 BLAST 的 ID 列表均生成后，THE ID_Merger SHALL 合并两个列表、排序去重，输出到 `Result_Dir/{Proj_Name}.identify.candidates.idlist`。
10. WHEN 合并 ID 列表生成后，THE ID_Merger SHALL 将 Transcript_ID 转换为 Gene_ID，输出到 `Result_Dir/{Proj_Name}.identify.candidates.gene.idlist`。
11. WHEN ID 列表生成后，THE Sequence_Extractor SHALL 从 CDS 数据库提取对应序列，输出到 `Result_Dir/{Proj_Name}.identify.candidates.cds.fa`。
12. WHEN ID 列表生成后，THE Sequence_Extractor SHALL 从蛋白数据库提取对应序列，输出到 `Result_Dir/{Proj_Name}.identify.candidates.pep.fa`。
13. WHEN Gene_ID 列表生成后，THE Sequence_Extractor SHALL 从基因组序列数据库提取对应序列，输出到 `Result_Dir/{Proj_Name}.identify.candidates.genome.fa`。
14. IF 参考蛋白文件和 HMM 文件均不存在，THEN THE Pipeline SHALL 输出错误信息提示用户提供至少一种输入，并以非零状态码退出。

---

### 需求 4：阶段二——系统发育树构建（tree）

**用户故事：** 作为研究人员，我希望基于候选成员蛋白序列自动构建 ML 系统发育树，以便分析基因家族的进化关系。

#### 验收标准

1. WHEN tree 阶段启动时，THE Tree_Builder SHALL 检查 `Result_Dir/{Proj_Name}.identify.candidates.pep.fa` 是否存在；IF 该文件不存在，THEN THE Pipeline SHALL 输出错误信息提示先运行 identify 阶段，并以非零状态码退出。
2. WHEN 蛋白序列文件存在时，THE Tree_Builder SHALL 调用 muscle 进行多序列比对，输出到 `Result_Dir/{Proj_Name}.tree.pep.afa`。
3. WHEN 多序列比对完成后，THE Tree_Builder SHALL 调用 trimal（`-automated1` 参数）修剪比对，输出到 `Result_Dir/{Proj_Name}.tree.pep.trimed.afa`。
4. WHEN 修剪后的比对文件生成后，THE Tree_Builder SHALL 调用 IQ-TREE（参数 `-bb {tree.iqtree_bootstrap} -bnni -nt AUTO`）构建 ML 树，输出文件前缀为 `Result_Dir/{Proj_Name}.tree.pep.trimed.afa`。
5. IF trimal 或 IQ-TREE 执行失败，THEN THE Pipeline SHALL 输出包含工具名称和错误信息的日志，并以非零状态码退出。

---

### 需求 5：阶段三——保守结构域分析（domain）

**用户故事：** 作为研究人员，我希望自动提交候选成员蛋白序列到 NCBI Batch CD-Search，以便批量获取保守结构域注释。

#### 验收标准

1. WHEN domain 阶段启动时，THE Domain_Analyzer SHALL 检查 `Result_Dir/{Proj_Name}.identify.candidates.pep.fa` 是否存在；IF 该文件不存在，THEN THE Pipeline SHALL 输出错误信息提示先运行 identify 阶段，并以非零状态码退出。
2. WHEN 蛋白序列文件存在时，THE Domain_Analyzer SHALL 向 NCBI Batch CD-Search API 提交序列，使用配置中的 `domain.evalue` 和 `domain.maxhit` 参数。
3. WHEN 搜索任务提交成功后，THE Domain_Analyzer SHALL 每隔 5 秒轮询任务状态，直到状态码为 0（完成）。
4. WHEN 搜索完成后，THE Domain_Analyzer SHALL 下载结果并保存到 `Result_Dir/{Proj_Name}.domain.cdd.txt`。
5. IF NCBI API 返回非 200 状态码，THEN THE Domain_Analyzer SHALL 输出包含 HTTP 状态码和原因的错误信息，并以非零状态码退出。
6. IF 任务状态码为 1、2、4 或 5（错误状态），THEN THE Domain_Analyzer SHALL 输出对应的错误描述，并以非零状态码退出。
7. WHILE 轮询任务状态时，THE Domain_Analyzer SHALL 在每次轮询时输出进度提示（如 "等待结果中..."）。

---

### 需求 6：阶段三（子阶段）——基于结构域筛选基因（domain-filter）

**用户故事：** 作为研究人员，我希望从全基因组蛋白中筛选出含有与目标基因家族相同或同一超家族结构域的基因，以便补充和验证候选成员集合。

#### 验收标准

1. WHEN domain-filter 阶段启动时，THE Domain_Filter SHALL 检查 `Result_Dir/{Proj_Name}.domain.cdd.txt` 是否存在；IF 该文件不存在，THEN THE Pipeline SHALL 输出错误信息提示先运行 domain 阶段，并以非零状态码退出。
2. WHEN `domain.target_domains` 在配置中已指定时，THE Domain_Filter SHALL 使用该列表作为目标结构域 accession；WHEN `domain.target_domains` 未指定时，THE Domain_Filter SHALL 自动从 `Result_Dir/{Proj_Name}.domain.cdd.txt` 中提取目标基因家族成员所含的 Domain_Accession 和 Superfamily_Accession 作为目标列表。
3. WHEN `domain.genome_cdd` 在配置中已指定且对应文件存在时，THE Domain_Filter SHALL 解析该全基因组 CDD 结果文件；IF `domain.genome_cdd` 未指定，THEN THE Domain_Filter SHALL 对全基因组蛋白库（`databases.pep`）分批提交 NCBI Batch CD-Search，每批不超过 4000 条序列，并将合并结果保存到 `Result_Dir/{Proj_Name}.domain.genome.cdd.txt`。
4. WHEN 全基因组 CDD 结果可用时，THE Domain_Filter SHALL 从中筛选含有目标 Domain_Accession（精确匹配）或目标 Superfamily_Accession（超家族匹配）的基因 ID，输出到 `Result_Dir/{Proj_Name}.domain-filter.candidates.idlist`。
5. WHEN 筛选完成后，THE Domain_Filter SHALL 生成包含结构域信息的汇总表，输出到 `Result_Dir/{Proj_Name}.domain-filter.summary.tsv`，列为：`gene_id`、`domain_accession`、`domain_name`、`superfamily_accession`、`evalue`。
6. IF 全基因组 CDD 结果文件解析失败，THEN THE Domain_Filter SHALL 输出包含文件路径和错误原因的错误信息，并以非零状态码退出。

---

### 需求 7：阶段四——Motif 发现与筛选（motif）

**用户故事：** 作为研究人员，我希望对候选基因家族成员运行 MEME 发现保守 motif，并用 FIMO 在全基因组蛋白中扫描，以便筛选出含有相同 motif 组合的基因。

#### 验收标准

1. WHEN motif 阶段启动时，THE Motif_Analyzer SHALL 检查 `Result_Dir/{Proj_Name}.identify.candidates.pep.fa` 是否存在；IF 该文件不存在，THEN THE Pipeline SHALL 输出错误信息提示先运行 identify 阶段，并以非零状态码退出。
2. WHEN 蛋白序列文件存在时，THE Motif_Analyzer SHALL 调用 meme，使用参数 `-nmotifs {motif.num_motifs} -minw {motif.min_width} -maxw {motif.max_width} -protein`，对候选成员蛋白序列进行 motif 发现，输出目录为 `Result_Dir/{Proj_Name}.motif.meme/`。
3. WHEN MEME 运行完成后，THE Motif_Analyzer SHALL 调用 fimo，使用 MEME 输出的 motif 文件扫描全基因组蛋白数据库（`databases.pep`），输出目录为 `Result_Dir/{Proj_Name}.motif.fimo/`。
4. WHEN FIMO 运行完成后，THE Motif_Analyzer SHALL 解析 FIMO 结果，统计每个基因命中的 motif 列表。
5. WHEN `motif.filter_mode` 为 `any` 时，THE Motif_Analyzer SHALL 将含有任意一个目标 motif 的基因纳入筛选结果。
6. WHEN `motif.filter_mode` 为 `all` 时，THE Motif_Analyzer SHALL 仅将含有全部目标 motif 的基因纳入筛选结果。
7. WHEN `motif.filter_mode` 为 `min_count` 时，THE Motif_Analyzer SHALL 仅将命中 motif 数量不少于 `motif.min_motif_count` 的基因纳入筛选结果。
8. WHEN 筛选完成后，THE Motif_Analyzer SHALL 将筛选后的基因 ID 列表输出到 `Result_Dir/{Proj_Name}.motif-filter.candidates.idlist`。
9. WHEN 筛选完成后，THE Motif_Analyzer SHALL 生成每个基因的 motif 命中汇总表，输出到 `Result_Dir/{Proj_Name}.motif-filter.summary.tsv`，列为：`gene_id`、`motif_id`、`motif_name`、`sequence_name`、`start`、`stop`、`score`、`p_value`。
10. IF meme 或 fimo 执行失败，THEN THE Pipeline SHALL 输出包含工具名称和错误信息的日志，并以非零状态码退出。

---

### 需求 8：阶段五——基因组共线性分析（collinearity）

**用户故事：** 作为研究人员，我希望分析基因家族成员在基因组中的共线性分布和复制模式，以便研究基因家族的扩张机制。

#### 验收标准

1. WHEN collinearity 阶段启动时，THE Collinearity_Analyzer SHALL 检查 `Result_Dir/{Proj_Name}.identify.candidates.gene.idlist`、`databases.gff3` 和 `databases.pep` 是否均存在；IF 任意文件不存在，THEN THE Pipeline SHALL 输出包含缺失文件路径的错误信息，并以非零状态码退出。
2. WHEN 所需输入文件均存在时，THE Collinearity_Analyzer SHALL 调用 blastp 对全基因组蛋白进行全对全自比对（all-vs-all），线程数为 `collinearity.blast_threads`，输出到 `Result_Dir/{Proj_Name}.collinearity.blast.out`。
3. WHEN `collinearity.tool` 为 `jcvi` 时，THE Collinearity_Analyzer SHALL 调用 JCVI（python -m jcvi.compara.catalog ortholog）进行共线性分析，输出共线性块文件到 `Result_Dir/{Proj_Name}.collinearity/`。
4. WHEN `collinearity.tool` 为 `mcscanx` 时，THE Collinearity_Analyzer SHALL 调用 MCScanX 进行共线性分析，输出共线性块文件到 `Result_Dir/{Proj_Name}.collinearity/`。
5. WHEN 共线性分析完成后，THE Collinearity_Analyzer SHALL 从共线性结果中提取包含目标基因家族成员的共线性块，输出到 `Result_Dir/{Proj_Name}.collinearity.blocks.tsv`，列为：`block_id`、`gene1`、`gene2`、`chromosome1`、`chromosome2`、`score`。
6. WHEN 共线性分析完成后，THE Collinearity_Analyzer SHALL 解析 GFF3 文件，生成基因家族成员在染色体上的位置信息表，输出到 `Result_Dir/{Proj_Name}.collinearity.gene-location.tsv`，列为：`gene_id`、`chromosome`、`start`、`end`、`strand`。
7. IF blastp 全对全自比对执行失败，THEN THE Pipeline SHALL 输出包含错误信息的日志，并以非零状态码退出。
8. IF 共线性分析工具执行失败，THEN THE Pipeline SHALL 输出包含工具名称和错误信息的日志，并以非零状态码退出。

---

### 需求 9：阶段六——基因属性分析（properties）

**用户故事：** 作为研究人员，我希望 CLI 能够计算候选成员的理化性质并提示在线分析步骤，以便我知道下一步该做什么。

#### 验收标准

1. WHEN properties 阶段启动时，THE Pipeline SHALL 检查 `Result_Dir/{Proj_Name}.identify.candidates.pep.fa` 是否存在；IF 该文件不存在，THEN THE Pipeline SHALL 输出错误信息提示先运行 identify 阶段，并以非零状态码退出。
2. WHEN properties 阶段执行时，THE Pipeline SHALL 使用 BioPython 计算候选成员蛋白的理化性质（分子量、等电点、氨基酸组成），并将结果输出到 `Result_Dir/{Proj_Name}.properties.tsv`。
3. WHEN 理化性质计算完成后，THE Pipeline SHALL 输出提示信息，告知用户跨膜域预测（DeepTMHMM）、信号肽预测（SignalP-6.0）和亚细胞定位预测（WoLF PSORT）的在线工具地址及所需输入文件路径。

---

### 需求 10：阶段七——转录组分析（trans）

**用户故事：** 作为研究人员，我希望 CLI 能够基于用户提供的表达量数据生成热图和差异基因列表，以便快速可视化基因家族的表达模式。

#### 验收标准

1. WHEN trans 阶段启动时，THE Pipeline SHALL 检查配置中 `trans.expression_matrix` 字段指定的表达量矩阵文件是否存在；IF 该文件不存在，THEN THE Pipeline SHALL 输出错误信息提示用户提供表达量矩阵，并以非零状态码退出。
2. WHEN 表达量矩阵文件存在时，THE Pipeline SHALL 筛选出属于当前基因家族成员的行，生成子矩阵。
3. WHEN 子矩阵生成后，THE Pipeline SHALL 使用 seaborn 或 matplotlib 绘制表达量热图，输出到 `Result_Dir/{Proj_Name}.trans.expression-heatmap.pdf`。
4. WHERE 配置中 `trans.logfc_threshold` 和 `trans.pvalue_threshold` 字段已设置，THE Pipeline SHALL 根据阈值筛选差异表达基因，输出到 `Result_Dir/{Proj_Name}.trans.deg.tsv`。
5. WHERE 配置中 `trans.venn_groups` 字段已设置（包含两个或多个样本组），THE Pipeline SHALL 绘制 Venn 图，输出到 `Result_Dir/{Proj_Name}.trans.venn.pdf`。
6. IF 表达量矩阵中没有任何基因家族成员，THEN THE Pipeline SHALL 输出警告信息并跳过热图绘制步骤。

---

### 需求 11：目录与文件管理

**用户故事：** 作为研究人员，我希望 pipeline 自动管理输入输出目录和文件命名，以便保持结果的一致性和可追溯性。

#### 验收标准

1. THE Pipeline SHALL 在每个阶段开始前自动创建所需的输出目录（若已存在则跳过，不报错）。
2. THE Pipeline SHALL 按照 `results/{Proj_Name}_results/` 的规则命名结果目录，其中 `{Proj_Name}` 来自配置文件的 `project_name` 字段。
3. THE Pipeline SHALL 按照 `data/{Proj_Name}_data/` 的规则命名数据目录。
4. THE Pipeline SHALL 对所有输出文件统一使用 `{Proj_Name}.{stage}.{type}.{ext}` 的命名格式，其中 `stage` 为阶段名称，`type` 为文件类型描述，`ext` 为文件扩展名。
5. WHEN 某个中间文件已存在时，THE Pipeline SHALL 默认跳过重新生成该文件的步骤，并输出跳过提示；WHEN `--force` 标志被设置时，THE Pipeline SHALL 强制重新生成所有中间文件。
6. THE Pipeline SHALL 在 `Result_Dir/pipeline.log` 中记录每次运行的时间戳、执行的阶段、每个步骤的状态（成功/跳过/失败）及耗时。

---

### 需求 12：全流程汇总报告

**用户故事：** 作为研究人员，我希望在全流程运行结束后获得一份汇总报告，以便快速了解各阶段的产出情况和关键统计数字。

#### 验收标准

1. WHEN `gfpipeline run` 全流程执行完成后，THE Summary_Reporter SHALL 生成汇总报告文件 `Result_Dir/{Proj_Name}.summary.txt`。
2. THE Summary_Reporter SHALL 在汇总报告中列出各阶段产出的关键文件路径，包括：identify 阶段的候选基因 ID 列表、tree 阶段的系统发育树文件、domain 阶段的 CDD 结果文件、domain-filter 阶段的筛选 ID 列表、motif 阶段的 MEME 输出目录和筛选 ID 列表、collinearity 阶段的共线性块文件和基因位置文件。
3. THE Summary_Reporter SHALL 在汇总报告中列出各阶段的关键统计数字，包括：identify 阶段的候选基因数量、domain-filter 阶段筛选后的基因数量、motif 阶段发现的 motif 数量和筛选后的基因数量、collinearity 阶段识别的共线性块数量。
4. IF 某阶段未执行或执行失败，THEN THE Summary_Reporter SHALL 在汇总报告中对应位置标注"未执行"或"执行失败"，而非跳过该条目。
5. FOR ALL 成功执行的阶段，THE Summary_Reporter 生成的汇总报告中列出的文件路径 SHALL 均指向实际存在的文件（文件存在性属性）。

---

### 需求 13：工具依赖检查

**用户故事：** 作为研究人员，我希望 pipeline 在运行前检查所需外部工具是否可用，以便在分析中途失败前得到明确提示。

#### 验收标准

1. WHEN Pipeline 启动时，THE Pipeline SHALL 检查当前阶段所需的所有外部工具是否在系统 PATH 或配置文件指定路径中可执行，各阶段所需工具如下：identify 阶段：muscle、hmmbuild、hmmsearch、hmmemit、blastp；tree 阶段：muscle、trimal、iqtree2；domain 阶段：无本地工具；motif 阶段：meme、fimo；collinearity 阶段：blastp、mcscanx 或 jcvi；genome-db blast 阶段：makeblastdb；genome-db gff-qc 阶段：samtools（必需），minimap2（当 `gff_qc.transcript_fa` 已配置时必需），stringtie（当 `gff_qc.rna_bam` 已配置时必需）；genome-db gene-index 阶段：samtools；genome-db rep-index 阶段：无外部工具（纯 Python 实现）。
2. IF 任意所需工具不可执行，THEN THE Pipeline SHALL 输出包含缺失工具名称和建议安装方式的错误信息，并以非零状态码退出。
3. THE Pipeline SHALL 仅检查当前将要运行的阶段所需的工具，而非所有阶段的工具。
4. WHEN genome-db gff-qc 阶段启动且 `gff_qc.transcript_fa` 已配置时，THE Pipeline SHALL 在工具检查阶段验证 minimap2 是否可执行；IF minimap2 不可执行，THEN THE Pipeline SHALL 输出包含缺失工具名称和建议安装方式的错误信息，并以非零状态码退出。
5. WHEN genome-db gff-qc 阶段启动且 `gff_qc.rna_bam` 已配置时，THE Pipeline SHALL 在工具检查阶段验证 stringtie 是否可执行；IF stringtie 不可执行，THEN THE Pipeline SHALL 输出包含缺失工具名称和建议安装方式的错误信息，并以非零状态码退出。

---

### 需求 14：基因组 BLAST 数据库构建（genome-db blast）

**用户故事：** 作为研究人员，我希望通过 `gfpipeline genome-db blast` 命令自动为基因组蛋白和核酸序列构建 BLAST 数据库，以便后续 blastp/blastn 搜索可以直接使用。

#### 验收标准

1. WHEN 用户执行 `gfpipeline genome-db blast --config <config_file>` 时，THE BLAST_DB_Builder SHALL 读取配置中 `databases.pep` 和 `databases.genome` 字段，确认输入文件存在后开始建库流程。
2. WHEN `genome_db.build_prot` 为 `true` 且 `databases.pep` 文件存在时，THE BLAST_DB_Builder SHALL 调用 `makeblastdb -dbtype prot` 为蛋白序列文件构建 BLAST 数据库，数据库文件前缀由 `databases.blast_db` 指定；若 `databases.blast_db` 未配置，则前缀与输入文件路径相同。
3. WHEN `genome_db.build_nucl` 为 `true` 且 `databases.genome` 文件存在时，THE BLAST_DB_Builder SHALL 调用 `makeblastdb -dbtype nucl` 为核酸序列文件构建 BLAST 数据库，数据库文件前缀与输入文件路径相同。
4. WHEN 建库完成后，THE BLAST_DB_Builder SHALL 输出统计信息，包括：序列数量、数据库类型（prot 或 nucl）、生成的数据库文件路径列表（`.phr`/`.pin`/`.psq` 等）。
5. IF `genome_db.build_prot` 和 `genome_db.build_nucl` 均为 `false`，THEN THE BLAST_DB_Builder SHALL 输出警告信息提示用户至少启用一种建库类型，并以非零状态码退出。
6. IF makeblastdb 执行失败，THEN THE BLAST_DB_Builder SHALL 输出包含工具名称、输入文件路径和错误原因的错误信息，并以非零状态码退出。
7. WHEN 对应的 BLAST 数据库文件已存在且 `--force` 标志未设置时，THE BLAST_DB_Builder SHALL 跳过该数据库的构建并输出跳过提示。

---

### 需求 15：GFF 文件质量检查与修正（genome-db gff-qc）

**用户故事：** 作为研究人员，我希望通过 `gfpipeline genome-db gff-qc` 命令对 GFF3 注释文件进行全面质检并修正疑似不完整的基因模型，以便为后续分析提供高质量的基因注释。

#### 验收标准

1. WHEN 用户执行 `gfpipeline genome-db gff-qc --config <config_file>` 时，THE GFF_QC SHALL 检查 `databases.gff3` 和 `databases.genome` 文件是否存在；IF 任意文件不存在，THEN THE GFF_QC SHALL 输出包含缺失文件路径的错误信息，并以非零状态码退出。
2. WHEN 输入文件均存在时，THE GFF_QC SHALL 检查 GFF3 格式规范性，包括：ID/Parent 属性完整性、坐标合法性（start ≤ end，坐标在染色体范围内）、strand 一致性（同一基因的所有子特征 strand 相同）。
3. WHEN GFF3 格式检查完成后，THE GFF_QC SHALL 检查每个基因模型的完整性，判定条件包括：是否具有 start codon 注释、是否具有 stop codon 注释、是否具有完整 CDS 记录。
4. WHEN 基因模型完整性检查完成后，THE GFF_QC SHALL 将满足以下任一条件的基因标记为 Truncated_Gene：位于 scaffold/contig 末端 `gff_qc.edge_distance` bp 范围内、CDS 总长度小于 `gff_qc.min_cds_len`、缺少 start codon 或 stop codon 注释。
5. WHEN 质检完成后，THE GFF_QC SHALL 将质检结果输出到 `{gff_qc.output_dir}/gff_qc.report.tsv`，列为：`gene_id`、`issue_type`、`detail`、`is_truncated`。
6. WHEN 质检完成后，THE GFF_QC SHALL 将质检统计摘要输出到 `{gff_qc.output_dir}/gff_qc.summary.txt`，包含：检查基因总数、各 issue_type 的数量、Truncated_Gene 数量。
7. WHEN `gff_qc.rna_bam` 已配置且对应 BAM 文件存在时，THE GFF_QC SHALL 使用 StringTie 对 Truncated_Gene 所在基因组区域重新预测转录本结构，并将预测结果与原注释合并，用于修正 CDS 边界。
8. WHEN `gff_qc.transcript_fa` 已配置且对应文件存在时，THE GFF_QC SHALL 调用 minimap2 将转录本序列比对到基因组，提取 Truncated_Gene 区域的比对结果，辅助修正 CDS 边界。
9. WHEN `gff_qc.extra_gff` 已配置且对应文件存在时，THE GFF_QC SHALL 将外部 GFF3 中与 Truncated_Gene 对应区域的注释与原注释合并，优先采用证据支持更充分的基因模型。
10. WHEN 修正步骤完成后，THE GFF_QC SHALL 将修正后的完整 GFF3 输出到 `{gff_qc.output_dir}/gff_qc.fixed.gff3`，未被修正的基因保持原始注释不变。
11. WHEN 修正步骤完成后，THE GFF_QC SHALL 将修正日志输出到 `{gff_qc.output_dir}/gff_qc.fix.log`，记录每个被修正基因的 `gene_id`、修正前注释摘要、修正后注释摘要及所用证据类型。
12. IF GFF3 文件解析失败，THEN THE GFF_QC SHALL 输出包含文件路径和错误原因的错误信息，并以非零状态码退出。
13. FOR ALL 未被标记为 Truncated_Gene 的基因，THE GFF_QC 输出的 `gff_qc.fixed.gff3` 中对应记录 SHALL 与原始 GFF3 中的记录完全一致（未修改基因保持不变属性）。

---

### 需求 16：基因序列索引构建（genome-db gene-index）

**用户故事：** 作为研究人员，我希望通过 `gfpipeline genome-db gene-index` 命令为基因组、CDS 和蛋白序列建立快速查询索引，以便通过基因 ID 或转录本 ID 高效提取序列。

#### 验收标准

1. WHEN 用户执行 `gfpipeline genome-db gene-index --config <config_file>` 时，THE Gene_Index_Builder SHALL 检查 `databases.gff3`（或 `gff_qc.output_dir}/gff_qc.fixed.gff3` 若存在）、`databases.genome`、`databases.cds`、`databases.pep` 是否均存在；IF 任意文件不存在，THEN THE Gene_Index_Builder SHALL 输出包含缺失文件路径的错误信息，并以非零状态码退出。
2. WHEN 输入文件均存在时，THE Gene_Index_Builder SHALL 使用 samtools faidx 或 BioPython SeqIO.index 为基因组序列（按染色体/scaffold）、CDS 序列（按 transcript ID）、蛋白序列（按 transcript ID）分别建立序列索引，索引文件输出到 `genome_db.index_dir`。
3. WHEN 序列索引构建完成后，THE Gene_Index_Builder SHALL 解析 GFF3 文件，构建基因 ID 到转录本 ID 的映射表，输出到 `{genome_db.index_dir}/gene2transcript.tsv`，列为：`gene_id`、`transcript_id`、`is_representative`。
4. WHEN 映射表构建完成后，THE Gene_Index_Builder SHALL 构建转录本 ID 到基因组坐标的映射表，输出到 `{genome_db.index_dir}/transcript2location.tsv`，列为：`transcript_id`、`chromosome`、`start`、`end`、`strand`、`exon_count`、`cds_length`。
5. WHEN 用户执行 `gfpipeline genome-db query --id <gene_or_transcript_id> --type [cds|pep|genome] --config <config_file>` 时，THE Gene_Index_Builder SHALL 利用已建立的索引快速提取对应序列，输出到 stdout；WHEN `--output <file>` 参数被指定时，THE Gene_Index_Builder SHALL 将序列输出到指定文件。
6. IF 查询的 ID 在索引中不存在，THEN THE Gene_Index_Builder SHALL 输出包含查询 ID 的错误信息，并以非零状态码退出。
7. FOR ALL 在 GFF3 中存在的 transcript ID，THE Gene_Index_Builder 构建的 `gene2transcript.tsv` 中 SHALL 包含该 transcript ID 对应的记录（索引完整性属性）。

---

### 需求 17：代表基因索引构建（genome-db rep-index）

**用户故事：** 作为研究人员，我希望通过 `gfpipeline genome-db rep-index` 命令自动为每个基因选取代表转录本并提取对应序列，以便后续分析使用统一的代表序列集合。

#### 验收标准

1. WHEN 用户执行 `gfpipeline genome-db rep-index --config <config_file>` 时，THE Rep_Index_Builder SHALL 检查 `databases.gff3`、`databases.cds`、`databases.pep` 是否均存在；IF 任意文件不存在，THEN THE Rep_Index_Builder SHALL 输出包含缺失文件路径的错误信息，并以非零状态码退出。
2. WHEN `genome_db.rep_selection` 为 `longest_cds`（默认）时，THE Rep_Index_Builder SHALL 对每个基因从其所有转录本中选取 CDS 序列最长的转录本作为 Representative_Transcript；若多个转录本 CDS 长度相同，则选取 transcript ID 字典序最小的转录本。
3. WHEN `genome_db.rep_selection` 为 `longest_mrna` 时，THE Rep_Index_Builder SHALL 对每个基因从其所有转录本中选取 mRNA 总长度最长的转录本作为 Representative_Transcript；若多个转录本 mRNA 长度相同，则选取 transcript ID 字典序最小的转录本。
4. WHEN 代表转录本选取完成后，THE Rep_Index_Builder SHALL 提取所有代表转录本的 CDS 序列，输出到 `{genome_db.rep_index_dir}/{Genome_Name}.rep.cds.fa`。
5. WHEN 代表转录本 CDS 序列提取完成后，THE Rep_Index_Builder SHALL 将代表转录本 CDS 序列翻译为蛋白序列，输出到 `{genome_db.rep_index_dir}/{Genome_Name}.rep.pep.fa`。
6. WHEN 代表转录本选取完成后，THE Rep_Index_Builder SHALL 输出代表转录本 ID 列表到 `{genome_db.rep_index_dir}/{Genome_Name}.rep.idlist`，每行一个 transcript ID。
7. WHEN 代表转录本选取完成后，THE Rep_Index_Builder SHALL 输出基因到代表转录本的映射表到 `{genome_db.rep_index_dir}/{Genome_Name}.gene2rep.tsv`，列为：`gene_id`、`rep_transcript_id`、`cds_length`、`transcript_count`。
8. FOR ALL 在 GFF3 中存在的基因，THE Rep_Index_Builder 生成的 `gene2rep.tsv` 中 SHALL 包含且仅包含该基因对应的一条记录（每基因唯一代表转录本属性）。
9. FOR ALL 在 `gene2rep.tsv` 中的 `rep_transcript_id`，THE Rep_Index_Builder 生成的 `{Genome_Name}.rep.cds.fa` 中 SHALL 包含对应的序列记录（代表序列完整性属性）。
10. IF GFF3 中某基因的所有转录本均无 CDS 记录，THEN THE Rep_Index_Builder SHALL 在日志中记录该基因 ID 及警告信息，并跳过该基因，不将其纳入输出文件。
