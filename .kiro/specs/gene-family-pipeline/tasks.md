# 实现计划：gene-family-pipeline

## 概述

将现有 Jupyter Notebook 基因家族分析工作流重构为结构化 Python CLI 应用。按模块从底层到上层逐步实现：核心工具层 → 配置层 → 基因组数据库准备层 → 各分析阶段 → CLI 入口层 → 汇总报告。

## 任务

- [x] 1. 初始化项目结构与核心工具层
  - 在 `GeneFamily-kiro/` 下创建 `gfpipeline/` 包及所有子包（`cli/`、`config/`、`stages/`、`genome_db/`、`core/`、`report/`），每个目录添加 `__init__.py`
  - 创建 `setup.py` 或 `pyproject.toml`，注册 `gfpipeline` 命令行入口点
  - 在 `core/logger.py` 中实现日志模块：同时写入 stderr 和 `Result_Dir/pipeline.log`，支持 INFO/DEBUG 级别切换
  - 在 `core/runner.py` 中实现 `ToolRunner` 类：封装 `subprocess` 调用，支持 `dry_run`/`verbose`，失败时抛出 `ExternalToolError`
  - 在 `core/` 中定义异常层次：`PipelineError`、`ConfigError`、`ToolNotFoundError`、`StageInputError`、`ExternalToolError`、`ApiError`
  - _需求：1.7、1.8、1.9、11.6_

- [x] 2. 实现配置层
  - [x] 2.1 在 `config/schema.py` 中实现所有 Pydantic 配置类（`ToolsConfig`、`DatabasesConfig`、`IdentifyConfig`、`TreeConfig`、`DomainConfig`、`MotifConfig`、`CollinearityConfig`、`GffQcConfig`、`GenomeDbConfig`、`TransConfig`、`PipelineConfig`），含默认值
    - _需求：2.1、2.2、2.3、2.4_
  - [ ]* 2.2 为配置往返解析编写属性测试
    - **属性 1：配置往返解析**
    - **验证：需求 2.6**
  - [x] 2.3 在 `config/loader.py` 中实现 `load_config(path)` 函数：读取 YAML 文件，用 Pydantic 验证，捕获 `ValidationError` 转换为 `ConfigError` 并报告缺失字段名
    - _需求：2.1、2.5_
  - [ ]* 2.4 为缺失必填字段时解析失败编写属性测试
    - **属性 2：缺失必填字段时解析失败并报告字段名**
    - **验证：需求 2.5**

- [x] 3. 实现核心序列工具与文件管理
  - [x] 3.1 在 `core/file_manager.py` 中实现 `FileManager` 类：`result()`、`data()`、`ensure_dirs()`、`skip_if_exists()` 方法
    - _需求：11.1、11.2、11.3、11.4、11.5_
  - [x] 3.2 在 `core/sequence.py` 中实现 `extract_fasta`（BioPython SeqIO.index）、`transcript_to_gene_id`、`parse_hmm_idlist`、`parse_blast_idlist` 四个函数
    - _需求：3.5、3.8、3.9、3.10、3.11、3.12、3.13_
  - [ ]* 3.3 为 ID 列表合并为并集且无重复编写属性测试
    - **属性 3：ID 列表合并为并集且无重复**
    - **验证：需求 3.9**
  - [ ]* 3.4 为 Transcript_ID 到 Gene_ID 转换格式正确性编写属性测试
    - **属性 4：Transcript_ID 到 Gene_ID 转换的格式正确性**
    - **验证：需求 3.10**
  - [x] 3.5 在 `core/tool_checker.py` 中实现 `check_tools(stage, config)` 和 `is_executable(tool_path)` 函数，含 `STAGE_TOOLS` 映射表及动态添加 jcvi/mcscanx 逻辑
    - _需求：13.1、13.2、13.3、13.4、13.5_

- [x] 4. 检查点——确保所有测试通过
  - 确保所有测试通过，如有问题请告知用户。

- [x] 5. 实现 genome-db 模块
  - [x] 5.1 在 `genome_db/blast_db.py` 中实现 `BlastDbBuilder` 类：检查输入文件、调用 `makeblastdb` 构建蛋白/核酸数据库、输出统计信息、支持跳过已存在数据库
    - _需求：14.1、14.2、14.3、14.4、14.5、14.6、14.7_
  - [x] 5.2 在 `genome_db/gff_qc.py` 中实现 `GffQc` 类：`parse_gff3`、`check_format`、`check_completeness`、`mark_truncated`、`fix_truncated`、`run` 方法；输出 `report.tsv`、`summary.txt`、`fixed.gff3`、`fix.log`
    - _需求：15.1–15.12_
  - [ ]* 5.3 为 Truncated_Gene 标记的充分条件编写属性测试
    - **属性 5：Truncated_Gene 标记的充分条件**
    - **验证：需求 15.4**
  - [ ]* 5.4 为未被标记为 Truncated_Gene 的基因在输出中保持不变编写属性测试
    - **属性 6：未被标记为 Truncated_Gene 的基因在输出中保持不变**
    - **验证：需求 15.13**
  - [x] 5.5 在 `genome_db/gene_index.py` 中实现 `GeneIndexBuilder` 类：`build_fasta_index`（samtools faidx / BioPython fallback）、`build_gene2transcript`、`build_transcript2location`、`query`、`run` 方法；输出 `gene2transcript.tsv`、`transcript2location.tsv`
    - _需求：16.1–16.6_
  - [ ]* 5.6 为 gene2transcript 映射完整性编写属性测试
    - **属性 7：gene2transcript 映射的完整性**
    - **验证：需求 16.7**
  - [x] 5.7 在 `genome_db/rep_index.py` 中实现 `RepIndexBuilder` 类：`select_representative`（longest_cds/longest_mrna 策略，等长取字典序最小）、`run` 方法；输出 `rep.cds.fa`、`rep.pep.fa`、`rep.idlist`、`gene2rep.tsv`
    - _需求：17.1–17.7、17.10_
  - [ ]* 5.8 为每个基因有且仅有一条代表转录本记录编写属性测试
    - **属性 8：每个基因有且仅有一条代表转录本记录**
    - **验证：需求 17.8**
  - [ ]* 5.9 为代表转录本序列完整性编写属性测试
    - **属性 9：代表转录本序列完整性**
    - **验证：需求 17.9**

- [x] 6. 检查点——确保所有测试通过
  - 确保所有测试通过，如有问题请告知用户。

- [x] 7. 实现 identify 阶段
  - [x] 7.1 在 `stages/identify.py` 中实现 `IdentifyStage` 类：`build_hmm`（muscle → hmmbuild）、`hmm_search`（hmmsearch → parse_hmm_idlist）、`blast_search`（hmmemit → blastp → parse_blast_idlist）、`merge_ids`、`extract_sequences`、`run` 方法
    - 使用 `FileManager` 管理所有输出路径，使用 `ToolRunner` 调用外部工具
    - 支持跳过已存在中间文件（`skip_if_exists`）
    - _需求：3.1–3.14_
  - [ ]* 7.2 为 identify 阶段核心逻辑编写单元测试
    - 测试 `merge_ids` 的去重与排序
    - 测试缺少 ref.fa 和 hmm 文件时的错误处理
    - _需求：3.14_

- [x] 8. 实现 tree 阶段
  - [x] 8.1 在 `stages/tree.py` 中实现 `TreeStage` 类：检查 `candidates.pep.fa` 存在性，依次调用 muscle → trimal（`-automated1`）→ iqtree2（`-bb {bootstrap} -bnni -nt AUTO`），输出 `tree.pep.afa`、`tree.pep.trimed.afa` 及 IQ-TREE 树文件
    - _需求：4.1–4.5_

- [x] 9. 实现 domain 与 domain-filter 阶段
  - [x] 9.1 在 `stages/domain.py` 中实现 `DomainStage` 类：`submit_batch`（POST 到 NCBI Batch CD-Search API）、`poll_status`（每 5 秒轮询，处理错误状态码）、`download_result`、`run` 方法；输出 `domain.cdd.txt`
    - _需求：5.1–5.7_
  - [x] 9.2 在 `stages/domain_filter.py` 中实现 `DomainFilterStage` 类：`extract_target_domains`、`parse_cdd_result`、`filter_by_domains`（精确匹配 + 超家族匹配）、`submit_genome_batch`（分批 ≤4000 条）、`run` 方法；输出 `domain-filter.candidates.idlist`、`domain-filter.summary.tsv`
    - _需求：6.1–6.6_
  - [ ]* 9.3 为结构域筛选结果正确性编写属性测试
    - **属性 10：结构域筛选结果的正确性**
    - **验证：需求 6.4**

- [x] 10. 实现 motif 阶段
  - [x] 10.1 在 `stages/motif.py` 中实现 `MotifStage` 类：`run_meme`（调用 meme，输出到 `motif.meme/`）、`run_fimo`（调用 fimo 扫描全基因组蛋白，输出到 `motif.fimo/`）、`parse_fimo`（解析 fimo.tsv，构建 gene_id → [motif_id] 映射）、`filter_genes`（支持 any/all/min_count 三种模式）、`run` 方法；输出 `motif-filter.candidates.idlist`、`motif-filter.summary.tsv`
    - _需求：7.1–7.10_
  - [ ]* 10.2 为 Motif 筛选模式正确性编写属性测试
    - **属性 11：Motif 筛选模式的正确性**
    - **验证：需求 7.5、7.6、7.7**

- [x] 11. 实现 collinearity、properties、trans 阶段
  - [x] 11.1 在 `stages/collinearity.py` 中实现 `CollinearityStage` 类：检查输入文件、`run_all_vs_all_blast`、`run_jcvi`/`run_mcscanx`（按配置选择）、`extract_target_blocks`、`parse_gene_locations`、`run` 方法；输出 `collinearity.blast.out`、`collinearity/`、`collinearity.blocks.tsv`、`collinearity.gene-location.tsv`
    - _需求：8.1–8.8_
  - [x] 11.2 在 `stages/properties.py` 中实现 `PropertiesStage` 类：使用 BioPython `ProteinAnalysis` 计算分子量、等电点、氨基酸组成，输出 `properties.tsv`，并打印在线工具提示信息
    - _需求：9.1–9.3_
  - [x] 11.3 在 `stages/trans.py` 中实现 `TransStage` 类：`load_expression_matrix`、`filter_family_members`、`plot_heatmap`（seaborn clustermap）、`filter_deg`、`plot_venn`（matplotlib_venn）、`run` 方法；输出 `trans.expression-heatmap.pdf`、`trans.deg.tsv`、`trans.venn.pdf`
    - _需求：10.1–10.6_

- [x] 12. 实现汇总报告模块
  - [x] 12.1 在 `report/summary.py` 中实现 `Summary_Reporter` 类：汇总各阶段关键文件路径和统计数字，对未执行或失败的阶段标注"未执行"/"执行失败"，输出 `{Proj}.summary.txt`
    - _需求：12.1–12.4_
  - [ ]* 12.2 为汇总报告中文件路径均实际存在编写属性测试
    - **属性 12：汇总报告中的文件路径均实际存在**
    - **验证：需求 12.5**

- [x] 13. 实现 CLI 入口层
  - [x] 13.1 在 `cli/main.py` 中用 click 定义顶层 `gfpipeline` 命令组，添加全局选项 `--config`/`-c`、`--dry-run`、`--verbose`、`--force`；在顶层捕获所有 `PipelineError`，输出格式化错误信息后以非零状态码退出
    - _需求：1.1、1.6、1.7、1.8、1.9_
  - [x] 13.2 在 `cli/cmd_genome_db.py` 中实现 `genome-db` 子命令组及 `blast`、`gff-qc`、`gene-index`、`rep-index`、`query` 子子命令；当未指定子子命令时按顺序执行全部四个子步骤
    - _需求：1.2、1.3_
  - [x] 13.3 在 `cli/cmd_run.py` 中实现 `run` 子命令：按顺序执行 identify → tree → domain → domain-filter → motif → collinearity，当 BLAST 数据库不存在时自动先执行 `genome-db blast`，全流程结束后调用 `Summary_Reporter`
    - _需求：1.4_
  - [x] 13.4 在 `cli/cmd_stages.py` 中实现各阶段独立子命令（`identify`、`tree`、`domain`、`domain-filter`、`motif`、`collinearity`、`properties`、`trans`），每个子命令在启动前调用 `check_tools`
    - _需求：1.5、13.1–13.5_

- [x] 14. 最终检查点——确保所有测试通过
  - 确保所有测试通过，如有问题请告知用户。

## 备注

- 标有 `*` 的子任务为可选项，可在快速 MVP 阶段跳过
- 每个任务均引用具体需求条款以保证可追溯性
- 检查点确保增量验证，及早发现问题
- 属性测试使用 `hypothesis` 库，每个属性为独立子任务
- 单元测试验证具体示例和边界条件，属性测试验证普遍性质
