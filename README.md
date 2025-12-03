# Using-PHLASH-Population-History-Analysis-Pipeline
a developed, robust, expandable pipeline to use phlash to analyse Ne curve and calculate Split time base on CCR method.
- Start from VCF files
- For Slurm environment on HPC

---


# **VCF_Extraction.sh** 
是一个多群体 VCF 预处理脚本，用于：

* 批量提取指定群体的样本
* 过滤为 biallelic SNP
* 进行次要等位基因计数 (MAC ≥ 2) 过滤
* 剔除 mask 区域（如低质量或难比对区域）
* 输出最终标准化 VCF，用于群体遗传分析（如 PHLASH、IBD、PCA 等）

特点：**通用、可扩展、可直接在 SLURM 作业环境中运行**。

---

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
