# gfpipeline — 基因家族分析流程工具

`gfpipeline` 是一个命令行驱动的基因家族分析流程，将 HMM 搜索、BLAST 比对、系统发育树构建、保守结构域分析、Motif 发现、共线性分析等步骤整合为一个可配置、可重复的自动化流程。

---

## 目录

- [安装](#安装)
- [快速开始](#快速开始)
- [配置文件详解](#配置文件详解)
- [命令参考](#命令参考)
  - [全局选项](#全局选项)
  - [run — 执行全流程](#run--执行全流程)
  - [genome-db — 基因组数据库准备](#genome-db--基因组数据库准备)
  - [identify — 基因家族鉴定](#identify--基因家族鉴定)
  - [tree — 系统发育树构建](#tree--系统发育树构建)
  - [domain — 保守结构域分析](#domain--保守结构域分析)
  - [domain-filter — 结构域筛选](#domain-filter--结构域筛选)
  - [motif — Motif 发现与筛选](#motif--motif-发现与筛选)
  - [motif-filter — 仅重跑 Motif 筛选](#motif-filter--仅重跑-motif-筛选)
  - [refine — 交互式二次筛选与分析](#refine--交互式二次筛选与分析)
  - [collinearity — 共线性分析](#collinearity--共线性分析)
  - [properties — 理化性质分析](#properties--理化性质分析)
  - [trans — 转录组分析](#trans--转录组分析)
- [输出文件说明](#输出文件说明)
- [外部工具依赖](#外部工具依赖)
- [常见问题](#常见问题)

---

## 安装

### 环境要求

- Python ≥ 3.10
- pip

### 安装步骤

```bash
# 克隆仓库
git clone <repo_url>
cd GFpipeline

# 安装（推荐使用虚拟环境）
pip install -e .
```

安装完成后，`gfpipeline` 命令即可在终端使用：

```bash
gfpipeline --help
```


---

## 快速开始

### 第一步：准备数据

在项目根目录下准备以下文件：

| 文件 | 说明 |
|------|------|
| `data/{GeneFamily}_data/{GeneFamily}.ref.fa` | 目标基因家族的参考蛋白序列（FASTA 格式） |
| `data/Oryza_sativa.dna_sm.toplevel.fa` | 全基因组序列 |
| `data/Oryza_sativa.58.chr.gff3` | 基因组 GFF3 注释文件 |

> **命名约定：** 参考蛋白文件建议命名为 `{GeneFamily}.ref.fa`（如 `wrky.ref.fa`、`ARF.ref.fa`）。
> 流程会自动从文件名中提取基因家族名（stem），并以此命名衍生文件（HMM profile、hmmemit 输出等）。
> 例如 `wrky.ref.fa` → `wrky.hmm`、`wrky.hmmemit.out`。

代表蛋白序列、代表 CDS 序列及 BLAST 数据库均由 `genome-db` 阶段自动生成，无需手动准备。

### 第二步：创建配置文件

在项目根目录创建 `config.yaml`（详见[配置文件详解](#配置文件详解)）。

### 第三步：执行全流程分析

```bash
gfpipeline --config config.yaml run
```

全流程按顺序执行：identify → tree → domain → domain-filter → motif → collinearity。

> 若 BLAST 数据库不存在，会在 identify 之前自动先执行 `genome-db` 建库，无需手动操作。

### 第四步（可选）：交互式二次筛选

初步分析完成后，根据进化树拓扑、结构域分布、motif 组合等结果，手动判断哪些基因是真正的家族成员，将精炼后的基因 ID 写入一个 idlist 文件：

```
# data/my_refined_genes.idlist
Os01g0001
Os02g0002
Os03t0003-01   # 也可以直接写转录本 ID
```

在配置文件中指定该文件路径，然后执行二次分析：

```bash
# config.yaml 中添加：
# refinement:
#   idlist: data/my_refined_genes.idlist

gfpipeline --config config.yaml refine
```

二次分析会对精炼后的序列重新执行 tree、domain、motif 三项分析，输出文件以 `{Proj}.refine.*` 命名。

### 第五步：查看结果

所有结果输出到配置文件中 `result_dir` 指定的目录（默认 `results/`），汇总报告见 `results/{GeneFamily}.summary.txt`。

---

## 配置文件详解

配置文件为 YAML 格式，默认文件名为 `config.yaml`。以下是完整示例及各字段说明。

```yaml
# ============================================================
# 必填字段
# ============================================================

# 项目名称，用于命名数据目录、结果目录及所有中间文件
project_name: MyGeneFamily

# 数据目录（存放参考序列、HMM 文件等）
data_dir: data/MyGeneFamily_data/

# 结果输出目录
result_dir: results/MyGeneFamily_results/

databases:
  ref_fa:  data/MyGeneFamily_data/wrky.ref.fa                  # 目标基因家族参考蛋白序列
  genome:  data/Oryza_sativa.IRGSP-1.0.dna_sm.toplevel.fa     # 全基因组序列 FASTA
  gff3:    data/Oryza_sativa.IRGSP-1.0.58.chr.gff3            # GFF3 注释文件

# ref_fa 文件名的 stem（去掉 .ref.fa 后缀）决定 HMM profile 等衍生文件的命名：
# wrky.ref.fa → wrky.hmm、wrky.hmmemit.out
# 代表蛋白序列、代表 CDS 序列及 BLAST 数据库由 genome-db 阶段自动生成，
# 默认输出到 data/rep_index/ 目录，无需手动配置。

# ============================================================
# 工具路径（可选，默认使用系统 PATH 中的工具）
# ============================================================
tools:
  muscle:      muscle       # 多序列比对
  hmmbuild:    hmmbuild     # HMM profile 构建
  hmmsearch:   hmmsearch    # HMM 搜索
  hmmemit:     hmmemit      # HMM consensus 序列生成
  blastp:      blastp       # 蛋白 BLAST 搜索
  trimal:      trimal       # 比对修剪
  iqtree:      iqtree2      # 系统发育树构建
  meme:        meme         # Motif 发现
  fimo:        fimo         # Motif 扫描
  mcscanx:     MCScanX      # 共线性分析（mcscanx 模式）
  makeblastdb: makeblastdb  # BLAST 数据库构建
  samtools:    samtools     # 序列索引
  minimap2:    minimap2     # 转录本比对（gff-qc 可选）
  stringtie:   stringtie    # 转录本重预测（gff-qc 可选）

# ============================================================
# 各阶段参数（可选，均有默认值）
# ============================================================

identify:
  hmm_evalue:    100     # HMM 搜索 e-value 阈值（宽松，用于初筛）
  blast_evalue:  1e-5    # BLAST 搜索 e-value 阈值
  blast_threads: 10      # BLAST 并行线程数

tree:
  iqtree_bootstrap: 1000  # IQ-TREE bootstrap 重复次数

domain:
  evalue:  0.01    # CD-Search e-value 阈值
  maxhit:  250     # 每条序列最大命中数
  # cdd_db: /path/to/cdd_db/Cdd              # 可选：本地 CDD 数据库路径（用于 rpsblast 本地模式）
  # genome_cdd: path/to/genome.cdd.txt       # 可选：已有的全基因组 CDD 结果，跳过重新提交
  # target_domains: "pfam00161 OR cl08249"   # 可选：目标结构域 accession（支持 AND/OR 逻辑）

motif:
  num_motifs:      10     # MEME 发现的 motif 数量
  min_width:       6      # motif 最小宽度（氨基酸数）
  max_width:       50     # motif 最大宽度（氨基酸数）
  fimo_pvalue:     1e-5   # FIMO p-value 阈值（越小越严格）
  filter_mode:     any    # 筛选模式：any（任意一个）| all（全部）| min_count（最少数量）
  min_motif_count: 1      # filter_mode=min_count 时的最少命中数

collinearity:
  tool:          jcvi   # 共线性工具：jcvi | mcscanx
  blast_threads: 10     # 全对全 BLAST 线程数

genome_db:
  build_prot:    true           # 是否构建蛋白 BLAST 数据库
  build_nucl:    true           # 是否构建核酸 BLAST 数据库
  index_dir:     data/genome_index/   # 基因序列索引输出目录
  rep_index_dir: data/rep_index/      # 代表转录本索引输出目录
  genome_name:   ~              # 基因组标识符（默认使用 project_name）
  rep_selection: longest_cds    # 代表转录本选取策略：longest_cds | longest_mrna

gff_qc:
  output_dir:    data/gff_qc/   # GFF 质检结果输出目录
  min_cds_len:   150            # CDS 最小长度（bp），低于此值标记为疑似不完整
  edge_distance: 500            # scaffold 末端距离阈值（bp）
  # rna_bam:                    # 可选：RNA-seq BAM 文件列表，用于修正不完整基因
  #   - data/sample1.bam
  # transcript_fa: ~            # 可选：转录本 FASTA，用于 minimap2 比对辅助修正
  # extra_gff: []               # 可选：外部 GFF3 文件列表，用于补充注释

trans:
  expression_matrix: ~          # 表达量矩阵文件路径（行为基因，列为样本）
  logfc_threshold:   ~          # 差异表达筛选 log2FC 阈值
  pvalue_threshold:  ~          # 差异表达筛选 p-value 阈值
  venn_groups: []               # Venn 图样本组列表

# ============================================================
# 交互式二次筛选（可选，初步分析完成后使用）
# ============================================================
# refinement:
#   idlist: data/my_refined_genes.idlist  # 用户手动筛选后的基因 ID 列表文件路径
```


---

## 命令参考

### 全局选项

所有子命令均支持以下全局选项：

| 选项 | 说明 |
|------|------|
| `--config` / `-c` | 配置文件路径（默认：当前目录下的 `config.yaml`） |
| `--dry-run` | 打印将要执行的命令，但不实际运行 |
| `--verbose` | 输出详细日志（包含完整命令行和文件路径） |
| `--force` | 强制重新生成所有中间文件（默认跳过已存在的文件） |

---

### run — 执行全流程

按顺序执行 identify → tree → domain → domain-filter → motif → collinearity 六个阶段。

```bash
gfpipeline --config config.yaml run
```

若 BLAST 数据库不存在，会在 identify 之前自动先执行 `genome-db blast`。

**适用场景：** 首次分析一个新的基因家族，一键完成所有步骤。

---

### genome-db — 基因组数据库准备

为后续分析准备标准化的基因组数据库。不指定子命令时，依次执行全部四个子步骤。

```bash
# 执行全部四个子步骤
gfpipeline --config config.yaml genome-db

# 仅构建 BLAST 数据库
gfpipeline --config config.yaml genome-db blast

# 仅执行 GFF3 质检
gfpipeline --config config.yaml genome-db gff-qc

# 仅构建基因序列索引
gfpipeline --config config.yaml genome-db gene-index

# 仅构建代表转录本索引
gfpipeline --config config.yaml genome-db rep-index

# 按 ID 查询序列
gfpipeline --config config.yaml genome-db query --id Os01g0936800 --type pep
gfpipeline --config config.yaml genome-db query --id Os01t0936800-01 --type cds --output out.fa
```

**genome-db blast** 子步骤说明：
- 读取 `databases.pep` 和 `databases.genome`，调用 `makeblastdb` 分别构建蛋白和核酸 BLAST 数据库
- 数据库前缀由 `databases.blast_db` 指定
- 若数据库文件已存在且未设置 `--force`，则跳过

**genome-db gff-qc** 子步骤说明：
- 对 GFF3 文件进行格式规范性检查（ID/Parent 完整性、坐标合法性、strand 一致性）
- 检查每个基因模型的完整性（start/stop codon、CDS 记录）
- 将疑似不完整基因（Truncated_Gene）标记并尝试修正
- 若配置了 `gff_qc.rna_bam`，使用 StringTie 辅助修正；若配置了 `gff_qc.transcript_fa`，使用 minimap2 辅助修正
- 输出质检报告、修正后的 GFF3 及修正日志

**genome-db gene-index** 子步骤说明：
- 使用 samtools faidx 为基因组、CDS、蛋白序列建立索引
- 构建基因 ID ↔ 转录本 ID 映射表（`gene2transcript.tsv`）
- 构建转录本 ID ↔ 基因组坐标映射表（`transcript2location.tsv`）

**genome-db rep-index** 子步骤说明：
- 为每个基因选取代表转录本（默认策略：CDS 最长；等长时取字典序最小的 transcript ID）
- 提取代表转录本的 CDS 和蛋白序列
- 输出代表 ID 列表和基因-代表转录本映射表

**genome-db query** 参数：

| 参数 | 说明 |
|------|------|
| `--id` | 基因 ID 或转录本 ID |
| `--type` | 序列类型：`cds` / `pep` / `genome` |
| `--output` | 输出文件路径（默认输出到 stdout） |

---

### identify — 基因家族鉴定

通过 HMM 搜索和 BLAST 搜索联合鉴定候选基因家族成员。

```bash
gfpipeline --config config.yaml identify
```

**前置条件：**
- `databases.ref_fa`（参考蛋白序列）或同目录下的 `{stem}.hmm`（已有 HMM profile）至少存在一个
  （`{stem}` 为 `ref_fa` 文件名去掉 `.ref.fa` 后缀，如 `wrky.ref.fa` → `wrky`）
- `genome-db` 阶段已完成（代表蛋白序列和 BLAST 数据库已生成）

**执行流程：**
1. 若 `.hmm` 不存在：muscle 多序列比对参考蛋白 → hmmbuild 构建 HMM profile
2. hmmsearch 在全基因组蛋白中搜索，提取 e-value < `identify.hmm_evalue` 的候选 ID
3. hmmemit 生成 consensus 序列 → blastp 搜索，提取候选 ID
4. 合并 HMM 和 BLAST 结果，去重
5. 将 transcript ID 转换为 gene ID
6. 提取候选成员的 CDS、蛋白、基因组序列

---

### tree — 系统发育树构建

基于候选成员蛋白序列构建 ML 系统发育树。

```bash
gfpipeline --config config.yaml tree
```

**前置条件：** identify 阶段已完成（`{Proj}.identify.candidates.pep.fa` 存在）

**执行流程：**
1. muscle 多序列比对候选蛋白序列
2. trimal（`-automated1`）修剪比对
3. IQ-TREE（`-bb {bootstrap} -bnni -nt AUTO`）构建 ML 树

---

### domain — 保守结构域分析

向 NCBI Batch CD-Search 提交候选蛋白序列，获取保守结构域注释。

```bash
gfpipeline --config config.yaml domain
```

**前置条件：** identify 阶段已完成

**说明：** 此步骤需要网络连接访问 NCBI API。提交后每 5 秒轮询一次任务状态，直到完成。

---

### domain-filter — 结构域筛选

从全基因组蛋白中筛选含有目标结构域的基因，补充和验证候选成员集合。

```bash
gfpipeline --config config.yaml domain-filter
```

**前置条件：** domain 阶段已完成（`{Proj}.domain.cdd.txt` 存在）

**筛选逻辑：**
- 若配置了 `domain.target_domains`，使用指定的结构域 accession 列表
- 否则自动从候选基因的 CDD 结果中提取目标结构域和超家族 accession
- 若配置了 `domain.genome_cdd`，直接使用已有的全基因组 CDD 结果；否则分批提交全基因组蛋白（每批 ≤ 4000 条）

---

### motif — Motif 发现与筛选

用 MEME 发现候选成员的保守 motif，再用 FIMO 在全基因组中扫描，筛选含相同 motif 的基因。

```bash
gfpipeline --config config.yaml motif
```

**前置条件：** identify 阶段已完成

若只想修改筛选参数（`filter_mode` / `min_motif_count`）而不重跑 MEME/FIMO，使用：

```bash
gfpipeline --config config.yaml motif-filter
```

**筛选模式（`motif.filter_mode`）：**

| 模式 | 说明 |
|------|------|
| `any` | 命中任意一个 motif 即纳入结果（默认） |
| `all` | 必须命中全部 motif 才纳入结果 |
| `min_count` | 命中 motif 数量 ≥ `motif.min_motif_count` 才纳入结果 |

---

### refine — 交互式二次筛选与分析

在初步分析（tree/domain/motif）完成后，根据结果手动筛选基因，对精炼后的序列集合重新执行 tree、domain、motif 三项分析。

**使用流程：**

1. 将筛选后的基因 ID 写入一个文本文件（每行一个 ID，支持基因 ID 或转录本 ID，`#` 开头为注释行）
2. 在配置文件中添加 `refinement.idlist` 字段指向该文件
3. 执行 `gfpipeline refine`

```bash
# 执行完整 refine 流程：序列提取 → tree → domain → motif
gfpipeline --config config.yaml refine

# 仅执行二次进化树分析
gfpipeline --config config.yaml refine-tree

# 仅执行二次结构域分析
gfpipeline --config config.yaml refine-domain

# 仅执行二次 motif 分析
gfpipeline --config config.yaml refine-motif
```

**前置条件：**
- identify、tree、domain、motif 阶段已完成（用于参考初步结果）
- 配置文件中已设置 `refinement.idlist` 字段
- idlist 文件已放置到指定路径

**说明：**
- 支持基因 ID（如 `Os01g0001`）和转录本 ID（如 `Os01t0001-01`）混合输入；基因 ID 会自动转换为对应的代表转录本 ID
- 未在数据库中找到的 ID 会记录警告日志，不中断流程
- 若只想重跑某个子阶段，可直接使用 `refine-tree`、`refine-domain`、`refine-motif`，但需确保 `{Proj}.refine.candidates.pep.fa` 已存在（即先运行过 `refine`）

---

### collinearity — 共线性分析

分析基因家族成员在基因组中的共线性分布和复制模式。

```bash
gfpipeline --config config.yaml collinearity
```

**前置条件：** identify 阶段已完成，且 `databases.gff3` 和 `databases.pep` 存在

**共线性工具（`collinearity.tool`）：**

| 工具 | 说明 |
|------|------|
| `jcvi`（默认） | 使用 JCVI（`python -m jcvi.compara.catalog ortholog`） |
| `mcscanx` | 使用 MCScanX |

---

### properties — 理化性质分析

计算候选成员蛋白的理化性质，并提示在线分析工具。

```bash
gfpipeline --config config.yaml properties
```

**前置条件：** identify 阶段已完成

**计算内容（纯本地，无需网络）：**

| 列名 | 说明 |
|------|------|
| `length` | 氨基酸序列长度 |
| `molecular_weight` | 分子量（Da） |
| `isoelectric_point` | 等电点（pI） |
| `instability_index` | 不稳定性指数（Guruprasad et al. 1990；< 40 为稳定蛋白） |
| `aliphatic_index` | 脂肪族指数（Ikai 1980；= Ala + 2.9×Val + 3.9×(Ile+Leu)） |
| `gravy` | 总平均疏水性 GRAVY（Kyte & Doolittle 1982） |
| `aromaticity` | 芳香族氨基酸比例（Phe+Trp+Tyr） |
| `A`…`Y` | 20 种标准氨基酸组成（%） |

**在线工具提示：** 跨膜域预测（DeepTMHMM）、信号肽预测（SignalP-6.0）、亚细胞定位预测（WoLF PSORT）

---

### trans — 转录组分析

基于用户提供的表达量矩阵，生成热图和差异基因列表。

```bash
gfpipeline --config config.yaml trans
```

**前置条件：** identify 阶段已完成，且 `trans.expression_matrix` 文件存在

**表达量矩阵格式：** TSV 或 CSV 文件，行为基因 ID，列为样本名称。

**可选功能：**
- 设置 `trans.logfc_threshold` 和 `trans.pvalue_threshold` 可输出差异表达基因列表
- 设置 `trans.venn_groups`（两个或多个样本组）可绘制 Venn 图


---

## 输出文件说明

所有结果输出到配置文件中 `result_dir` 指定的目录（默认 `results/`），数据文件输出到 `data_dir` 指定的目录（默认 `data/`）。

### identify 阶段

| 文件 | 说明 |
|------|------|
| `{Proj}.identify.ref.afa` | 参考蛋白多序列比对结果 |
| `{Proj}.identify.hmm.out` | hmmsearch 原始输出 |
| `{Proj}.identify.hmm.idlist` | HMM 搜索候选 ID 列表 |
| `{Proj}.identify.blast.out` | blastp 原始输出 |
| `{Proj}.identify.blast.idlist` | BLAST 搜索候选 ID 列表 |
| `{Proj}.identify.candidates.idlist` | 合并后候选转录本 ID 列表 |
| `{Proj}.identify.candidates.gene.idlist` | 候选基因 ID 列表 |
| `{Proj}.identify.candidates.cds.fa` | 候选成员 CDS 序列 |
| `{Proj}.identify.candidates.pep.fa` | 候选成员蛋白序列 |
| `{Proj}.identify.candidates.genome.fa` | 候选成员基因组序列 |

### tree 阶段

| 文件 | 说明 |
|------|------|
| `{Proj}.tree.pep.afa` | 蛋白多序列比对结果 |
| `{Proj}.tree.pep.trimed.afa` | trimal 修剪后的比对 |
| `{Proj}.tree.pep.trimed.afa.treefile` | IQ-TREE 最优树（Newick 格式） |
| `{Proj}.tree.pep.trimed.afa.contree` | 共识树 |

### domain 阶段

| 文件 | 说明 |
|------|------|
| `{Proj}.domain.cdd.txt` | 候选成员 NCBI CD-Search 结果 |

### domain-filter 阶段

| 文件 | 说明 |
|------|------|
| `data/cdd/{Genome_Name}.genome.cdd.txt` | 全基因组 CD-Search 结果（若未提供现有结果） |
| `data/cdd/{Genome_Name}.genome.cdd.rpsblast.tsv` | rpsblast 原始输出（本地模式） |
| `{Proj}.domain-filter.candidates.idlist` | 结构域筛选后的候选基因 ID 列表 |
| `{Proj}.domain-filter.summary.tsv` | 结构域筛选汇总表（含结构域信息） |

### motif 阶段

| 文件 | 说明 |
|------|------|
| `{Proj}.motif.meme/` | MEME 输出目录（含 meme.xml、logo 图等） |
| `{Proj}.motif.fimo/` | FIMO 扫描输出目录（含 fimo.tsv） |
| `{Proj}.motif-filter.candidates.idlist` | Motif 筛选后的候选基因 ID 列表 |
| `{Proj}.motif-filter.summary.tsv` | 每个基因的 motif 命中汇总表 |

### refine 阶段

| 文件 | 说明 |
|------|------|
| `{Proj}.refine.candidates.pep.fa` | 从 Rep_Pep_DB 提取的精炼蛋白序列 |
| `{Proj}.refine.tree.pep.afa` | 精炼蛋白多序列比对结果 |
| `{Proj}.refine.tree.pep.trimed.afa` | trimal 修剪后的比对 |
| `{Proj}.refine.tree.pep.trimed.afa.treefile` | 二次 IQ-TREE 最优树（Newick 格式） |
| `{Proj}.refine.domain.cdd.txt` | 精炼成员 NCBI CD-Search 结果 |
| `{Proj}.refine.motif.meme/` | 二次 MEME 输出目录 |
| `{Proj}.refine.motif.fimo/` | 二次 FIMO 扫描输出目录 |
| `{Proj}.refine.motif-filter.candidates.idlist` | 二次 motif 筛选后的候选基因 ID 列表 |
| `{Proj}.refine.motif-filter.summary.tsv` | 二次 motif 命中汇总表 |

### collinearity 阶段

| 文件 | 说明 |
|------|------|
| `data/collinearity_mcscanx/{Genome_Name}.genome.blast.out` | 全基因组蛋白全对全 BLAST 结果 |
| `data/collinearity_mcscanx/{Genome_Name}.gff` | MCScanX 输入 GFF（全基因组所有基因坐标） |
| `data/collinearity_mcscanx/{Genome_Name}.collinearity` | MCScanX 全基因组共线性结果 |
| `data/collinearity_mcscanx/{Genome_Name}.tandem` | MCScanX 全基因组串联重复结果 |
| `data/collinearity_mcscanx/{Genome_Name}.html/` | MCScanX 可视化输出 |
| `{Proj}.collinearity.blocks.tsv` | 包含目标基因的共线性块 |
| `{Proj}.collinearity.gene-location.tsv` | 基因家族成员染色体位置信息 |

### properties 阶段

| 文件 | 说明 |
|------|------|
| `{Proj}.properties.tsv` | 候选成员理化性质表（长度、分子量、等电点、不稳定性指数、脂肪族指数、GRAVY、芳香族比例、氨基酸组成） |

### trans 阶段

| 文件 | 说明 |
|------|------|
| `{Proj}.trans.expression-heatmap.pdf` | 表达量热图 |
| `{Proj}.trans.deg.tsv` | 差异表达基因列表（需配置阈值） |
| `{Proj}.trans.venn.pdf` | Venn 图（需配置 venn_groups） |

### genome-db 阶段

| 文件 | 说明 |
|------|------|
| `data/genome_index/gene2transcript.tsv` | 基因 ID ↔ 转录本 ID 映射表 |
| `data/genome_index/transcript2location.tsv` | 转录本 ID ↔ 基因组坐标映射表 |
| `data/rep_index/{Genome_Name}.rep.cds.fa` | 代表转录本 CDS 序列 |
| `data/rep_index/{Genome_Name}.rep.pep.fa` | 代表转录本蛋白序列 |
| `data/rep_index/{Genome_Name}.rep.idlist` | 代表转录本 ID 列表 |
| `data/rep_index/{Genome_Name}.gene2rep.tsv` | 基因 ↔ 代表转录本映射表 |
| `data/gff_qc/gff_qc.report.tsv` | GFF3 质检详细报告 |
| `data/gff_qc/gff_qc.summary.txt` | GFF3 质检统计摘要 |
| `data/gff_qc/gff_qc.fixed.gff3` | 修正后的 GFF3 文件 |
| `data/gff_qc/gff_qc.fix.log` | 基因修正日志 |

### 全局文件

| 文件 | 说明 |
|------|------|
| `results/pipeline.log` | 每次运行的完整日志（时间戳、阶段状态、耗时） |
| `results/{Proj}.summary.txt` | 全流程汇总报告（各阶段产出文件路径和统计数字） |

---

## 外部工具依赖

各阶段所需的外部工具如下。工具需在系统 PATH 中可执行，或在配置文件 `tools` 字段中指定完整路径。

| 阶段 | 必需工具 | 可选工具 |
|------|----------|----------|
| identify | muscle, hmmbuild, hmmsearch, hmmemit, blastp | — |
| tree | muscle, trimal, iqtree2 | — |
| domain | — （仅需网络连接） | — |
| domain-filter | — （仅需网络连接） | — |
| motif | meme, fimo | — |
| collinearity | blastp, jcvi 或 MCScanX | — |
| properties | — （纯 Python） | — |
| trans | — （纯 Python） | — |
| genome-db blast | makeblastdb | — |
| genome-db gff-qc | samtools | minimap2（配置 transcript_fa 时），stringtie（配置 rna_bam 时） |
| genome-db gene-index | samtools | — |
| genome-db rep-index | — （纯 Python） | — |

Pipeline 在每个阶段启动时会自动检查所需工具是否可执行，若缺失会输出明确的错误信息和建议安装方式。

### 工具安装参考

```bash
# HMMER（hmmbuild, hmmsearch, hmmemit, hmmemit）
conda install -c bioconda hmmer

# BLAST+（blastp, makeblastdb）
conda install -c bioconda blast

# muscle（多序列比对）
conda install -c bioconda muscle

# trimal（比对修剪）
conda install -c bioconda trimal

# IQ-TREE 2（系统发育树）
conda install -c bioconda iqtree

# MEME Suite（meme, fimo）
conda install -c bioconda meme

# samtools
conda install -c bioconda samtools

# JCVI（共线性分析）
pip install jcvi

# MCScanX（共线性分析，可选）
# 从 https://github.com/wyp1125/MCScanX 下载编译
```

---

## 常见问题

**Q：运行时提示"No such option: --config"**

`--config` 是全局选项，必须放在子命令之前：
```bash
# 正确
gfpipeline --config config.yaml run
gfpipeline --config config.yaml run --dry-run

# 错误
gfpipeline run --config config.yaml
```

**Q：运行时提示"找不到配置文件"**

确认当前目录下存在 `config.yaml`，或使用 `--config` 参数指定路径：
```bash
gfpipeline --config /path/to/my_config.yaml run
```

**Q：如何跳过已完成的步骤重新运行某个阶段？**

默认情况下，若中间文件已存在则跳过。若需强制重新生成，使用 `--force`：
```bash
gfpipeline --config config.yaml --force identify
```

**Q：domain 阶段运行很慢**

domain 阶段需要将序列提交到 NCBI 服务器并等待结果，耗时取决于序列数量和服务器负载，通常需要几分钟到几十分钟。可以先运行其他不依赖 domain 结果的阶段（如 tree、motif）。

**Q：如何使用已有的全基因组 CDD 结果，跳过重新提交？**

在配置文件中设置 `domain.genome_cdd` 字段指向已有的结果文件：
```yaml
domain:
  genome_cdd: data/my_genome.cdd.txt
```

**Q：如何查看某次运行的详细日志？**

使用 `--verbose` 标志，或查看 `results/{Proj}_results/pipeline.log` 文件：
```bash
gfpipeline --config config.yaml --verbose run
```

**Q：trans 阶段的表达量矩阵格式是什么？**

TSV 或 CSV 文件，第一列为基因 ID，其余列为各样本的表达量值（如 TPM、FPKM 或 read count）。示例：

```
gene_id    sample1    sample2    sample3
Os01g0001  12.5       8.3        15.2
Os01g0002  0.0        2.1        0.5
```

**Q：如何只分析部分阶段？**

直接调用对应的子命令即可，例如只运行 motif 分析：
```bash
gfpipeline --config config.yaml motif
```

注意：各阶段有前置依赖，请确保前置阶段已完成（详见各阶段的"前置条件"说明）。

**Q：`--dry-run` 有什么用？**

`--dry-run` 会打印所有将要执行的命令，但不实际运行，适合在正式运行前检查命令是否正确：
```bash
gfpipeline --config config.yaml --dry-run run
```

**Q：如何进行交互式二次筛选？**

初步分析完成后，根据进化树、结构域、motif 结果手动筛选基因，将基因 ID 写入 idlist 文件，在配置文件中添加 `refinement.idlist` 字段，然后运行：
```bash
gfpipeline --config config.yaml refine
```
支持基因 ID 和转录本 ID 混合输入。若只需重跑某个子分析，可单独调用 `refine-tree`、`refine-domain` 或 `refine-motif`。
