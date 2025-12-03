# Using-PHLASH-Population-History-Analysis-Pipeline
a developed, robust, expandable pipeline to use phlash to analyse Ne curve and calculate Split time base on CCR method.

## plot_Ne_CCRSplit_20251203
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
