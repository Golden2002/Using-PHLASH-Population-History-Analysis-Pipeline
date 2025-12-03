# Using-PHLASH-Population-History-Analysis-Pipeline
a developed, robust, expandable pipeline to use phlash to analyse Ne curve and calculate Split time base on CCR method.
- Start from VCF files
- For Slurm environment on HPC

---

# General Notice
**Chromosome naming convention:**

The pipeline uses GRCh38 reference genome coordinates. Chromosomes should be named as follows:

* Autosomes: `chr1`, `chr2`, …, `chr22` (or `chr1_23` if you combine all autosomes)
* Sex chromosomes: `chrX`, `chrY`
* Mitochondrial genome: `chrM`

> Notice that all positions and windows in the scripts are based on the GRCh38 coordinate system.

---

# **VCF_Extraction.sh** 
是一个多群体 VCF 预处理脚本，用于：

* 批量提取指定群体的样本
* 过滤为 biallelic SNP
* 进行次要等位基因计数 (MAC ≥ 2) 过滤
* 剔除 mask 区域（如低质量或难比对区域）
* 输出最终标准化 VCF，用于群体遗传分析（如 PHLASH、IBD、PCA 等）

特点：**通用、可扩展、可直接在 SLURM 作业环境中运行**。



## 使用说明

### 1. 修改参数

在脚本开头修改：

```bash
POPULATIONS=("POP1" "POP2")       # 目标群体列表
VCF="/path/to/input.vcf.gz"       # 输入VCF文件
MASK_FILE="/path/to/mask.bed.gz"  # mask区域文件
INFO_FILE="/path/to/info.txt"     # 样本信息文件，格式: 样本ID 群体名
OUTBASE="./output"                # 输出目录
```

### 2. 提交作业

在 SLURM 环境中直接运行：

```bash
sbatch VCF_Extraction.sh
```

### 3. 输出文件

| 文件名                      | 说明                      |
| ------------------------ | ----------------------- |
| `ALL_TARGET_SAMPLES.txt` | 提取的所有目标样本ID             |
| `step1.filtered.vcf.gz`  | 样本和MAC过滤后的VCF           |
| `final_filtered.vcf.gz`  | 最终VCF（去除mask区域，可用于下游分析） |

---

# run_phlash_CCR_20251202.sh

本脚本用于在 **HPC 集群（SLURM）环境** 下批量执行 **PHLASH（Pairwise Haplotype Lengths for Approximate Skyline History）群体历史推断分析**。支持：

* ✔️ 单个群体分析
* ✔️ 两个群体的合并分析（如 *A::B*）
* ✔️ 自动抽样
* ✔️ 自动生成 Python 分析脚本
* ✔️ 自动创建输出目录

适用于需要对多个群体、或多个群体组合进行高效批处理的场景。


## 功能概述

* 根据 SLURM `--array` 设置，自动判断执行：

  * **单群体模式**（如 `Jino`）
  * **合并群体模式**（如 `Jino::Han_N`）
* 从 INFO 文件中自动抽取样本，并按群体或组合分配
* 自动构建 PHLASH 所需 contig
* 对 **chr1–chr22** 进行全染色体批量分析
* 分别输出：

  * 📁 **models/**：PHLASH 模型文件
  * 📁 **plots/**：拟合图
  * 📁 **logs/**：运行日志
  * 📁 **python_scripts/**：本次任务的 Python 脚本备份

---

## 输入文件要求

| 文件类型        | 描述                          |
| ----------- | --------------------------- |
| `INFO_FILE` | 含样本与群体信息的文件，要求至少包含：样本ID、群体名 |
| `VCF`       | 已 phased 的全基因组 VCF          |
| `MASK_FILE` | 屏蔽区域 BED 文件（StrictMask）     |
| `SEED`      | 用于一致性抽样的随机种子文件              |

示例（可自定义）：

```bash
INFO_FILE="/path/to/info.txt"
VCF="/path/to/phased.vcf.gz"
MASK_FILE="/path/to/mask.bed.gz"
SEED="/path/to/random_seed.txt"
```

---

## 主要可修改参数

你可以根据需求调整脚本顶部的参数：

```bash
POPULATIONS=("Jino" "Han_N" "Tibetan")   # 单群体
COMBINED_POPULATIONS=("Jino::Han_N")     # 合并群体
SAMPLE_PER_POP=10                        # 每群抽样数量
MUTATION_RATE="1.25e-8"                  # 突变率
KNOTS=20                                 # PHLASH 模型节点数
REG_PENALTY=6.0                          # 正则系数
```

---

## 使用方法

### 1. 提交单群体与合并群体任务

脚本使用 SLURM 数组任务管理运行模式：

```
#SBATCH --array=0-3
```

若定义了：

```bash
POPULATIONS=("A" "B" "C")
COMBINED_POPULATIONS=("A::B")
```

那么：

* `task_id=0` → 群体 A
* `task_id=1` → 群体 B
* `task_id=2` → 群体 C
* `task_id=3` → 合并群体 A::B

你可以直接提交：

```bash
sbatch run_phlash.sh
```

---

## 分析流程说明

脚本的执行逻辑如下：

1. **根据 SLURM_ARRAY_TASK_ID 判断分析模式**
2. **从 INFO 文件中自动抽取样本**

   * 单群体：抽取 N 个
   * 合并群体：各抽 N/2 个
3. **构建输出目录**
4. **自动生成对应的 Python PHLASH 脚本**
5. **Python 脚本内：**

   * 遍历 chr1–chr22
   * 按窗口（1Mb）构建 contig
   * 过滤变异数不足的片段
   * 在所有有效 contig 上运行 PHLASH
6. **保存模型与图形**

---

## 输出结构

提交任务后，每个群体或合并群体都会生成：

```
results/
  └── POP_NAME/
       ├── models/          # PHLASH 模型参数
       ├── plots/           # skyline 图和拟合曲线
       ├── logs/            # 运行日志
       └── python_scripts/  # 自动生成的 phlash_xxx.py
```

示例：

```
results/Jino/
results/Jino::Han_N/
```

---

## 注意事项

* 请确保输入 VCF 已 phased
* 合并群体分析时，样本量不足会自动跳过空群体
* 若首次测试，建议把 Python 脚本里的染色体列表缩小：

```python
chromosomes = ["chr1", "chr2"]
```

---


# run_phlash_20251127_XYmt.sh

### PHLASH 性染色体分析脚本

**功能**：
该脚本用于对性染色体（X、Y、MT）进行 PHLASH 分析，功能包括：

* 群体样本提取
* 男性样本选择（仅针对Y染色体）
* 随机选择样本子集
* 性染色体 contig 构建
* PHLASH 模型拟合和绘图

**特点**：

* 支持多群体分析
* 可扩展到任意群体和染色体
* 输出结构化目录：`plots/`、`models/`、`logs/`、`python_scripts/`

---

### 使用方法

1. **修改参数**
   在脚本开头修改以下路径和参数：

```bash
POPULATIONS=("POP1" "POP2")
VCF_X="/path/to/chrX.vcf.gz"
VCF_Y="/path/to/chrY.vcf.gz"
VCF_MT="/path/to/chrM.vcf.gz"
MALE_SAMPLES_FILE="/path/to/male_samples.txt"
INFO_FILE="/path/to/sample_info.txt"
SEED="/path/to/random_seed.txt"
```

2. **提交SLURM作业**

```bash
sbatch PHLASH_SexChrom.sh
```

3. **输出文件结构**

```
results/
└─ POP_CHROM/
   ├─ plots/                # PHLASH绘图
   ├─ models/               # 拟合结果
   ├─ logs/                 # 日志
   └─ python_scripts/       # 生成的Python脚本
```

4. **注意事项**

* Y染色体分析仅使用男性样本，默认最多抽取15个样本
* X和MT染色体分析使用随机样本子集
* Python 脚本需根据需求进一步调用 `phlash.contig` 和 `phlash.fit`

---




# plot_Ne_CCRSplit_20251203
**Python绘图可视化脚本使用示例**：
- `--base_dir`：指定包含群体数据的基础目录路径
- `--output_dir`：指定输出图表和结果文件的目录路径
- `--pop1` / `--pop2`：指定要分析的两个群体名称
- `--batch`：启用批量分析模式
- `--pop_pairs`：批量分析的群体对列表
- `--config`：使用配置文件指定所有参数
  
```bash
# 单对分析（使用默认目录）
python ccr_analysis.py --pop1 Jino --pop2 Han_N

# 单对分析（指定基础目录）
python ccr_analysis.py --pop1 Jino --pop2 Han_N --base_dir /path/to/phlash/results

# 单对分析（指定基础目录和输出目录）
python ccr_analysis.py --pop1 Jino --pop2 Han_N \
  --base_dir /path/to/phlash/results \
  --output_dir /path/to/output

# 批量分析（使用默认目录）
python ccr_analysis.py --batch --pop_pairs "Jino:Han_N" "Han_N:Tibetan"

# 批量分析（指定基础目录）
python ccr_analysis.py --batch --base_dir /path/to/phlash/results \
  --pop_pairs "Jino:Han_N" "Han_N:Tibetan"

# 批量分析（指定基础目录和输出目录）
python ccr_analysis.py --batch \
  --base_dir /path/to/phlash/results \
  --output_dir /path/to/output \
  --pop_pairs "Jino:Han_N" "Han_N:Tibetan"

# 使用配置文件（指定基础目录和输出目录）
python ccr_analysis.py --batch --config config.json

# 显示详细帮助信息
python ccr_analysis.py -h
```
