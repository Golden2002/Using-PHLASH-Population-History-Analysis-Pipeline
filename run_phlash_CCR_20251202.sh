#!/bin/bash
#SBATCH --job-name=PHLASH_Analysis
#SBATCH --output=results/log/PHLASH_%A_%a.log
#SBATCH --error=results/log/PHLASH_%A_%a.err.log
#SBATCH --ntasks=1
#SBATCH --partition=batch
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --array=3-3  # 索引：0-2为单个群体，3为合并分析

# ===== 配置加载 =====
if [ -f "config.sh" ]; then
    source config.sh
else
    echo "错误: 找不到配置文件 config.sh"
    echo "请创建 config.sh 并设置以下变量:"
    echo "  VCF, INFO_FILE, PHLASH_ENV"
    exit 1
fi

# ===== 环境设置 =====
if [ -n "$PHLASH_ENV" ]; then
    source "$PHLASH_ENV"
fi

# ===== 用户可修改参数 =====
# 群体列表 - 根据实际研究修改
POPULATIONS=("Population1" "Population2" "Population3")
COMBINED_POPULATIONS=("Population1::Population2")  # 合并群体定义

# 检查必要变量
if [ -z "$VCF" ]; then
    echo "错误: VCF路径未在config.sh中设置"
    exit 1
fi

if [ -z "$INFO_FILE" ]; then
    echo "错误: 样本信息文件未在config.sh中设置"
    exit 1
fi

# 默认参数
MUTATION_RATE=${MUTATION_RATE:-"1.25e-8"}
KNOTS=${KNOTS:-20}
REG_PENALTY=${REG_PENALTY:-6.0}
SPLINE_TYPE=${SPLINE_TYPE:-"piecewise"}
SEED=${SEED:-"shuf_random_seed.txt"}
SAMPLE_PER_POP=${SAMPLE_PER_POP:-10}
MIN_CONTIGS=${MIN_CONTIGS:-50}  # 最小contig数阈值
MAX_CHROMOSOMES=${MAX_CHROMOSOMES:-22}  # 最大分析的染色体数

# ===== 根据任务索引决定分析类型 =====
if [ $SLURM_ARRAY_TASK_ID -lt ${#POPULATIONS[@]} ]; then
    # 单个群体分析
    POP=${POPULATIONS[$SLURM_ARRAY_TASK_ID]}
    IS_COMBINED=false
    ANALYSIS_TYPE="single"
    echo "=== 单个群体分析 ==="
    echo "群体: $POP"
else
    # 合并群体分析
    COMBINED_IDX=$((SLURM_ARRAY_TASK_ID - ${#POPULATIONS[@]}))
    
    if [ $COMBINED_IDX -ge ${#COMBINED_POPULATIONS[@]} ]; then
        echo "警告: 任务ID $SLURM_ARRAY_TASK_ID 超出合并群体范围"
        exit 0
    fi
    
    COMBINED_POP=${COMBINED_POPULATIONS[$COMBINED_IDX]}
    IS_COMBINED=true
    ANALYSIS_TYPE="combined"
    
    echo "=== 合并群体分析 ==="
    echo "合并群体: $COMBINED_POP"
    
    # 解析合并群体名称
    IFS='::' read -ra POP_ARRAY <<< "$COMBINED_POP"
    POP1=${POP_ARRAY[0]}
    POP2=${POP_ARRAY[1]}
    echo "包含群体: $POP1 和 $POP2"
fi

# ===== 提取样本 =====
if [ "$IS_COMBINED" = false ]; then
    # 单个群体：抽取指定数量的样本
    if [ ! -f "$INFO_FILE" ]; then
        echo "错误: 样本信息文件不存在: $INFO_FILE"
        exit 1
    fi
    
    POP_SAMPLES=$(awk -v pop="$POP" '$2==pop{print $1}' "$INFO_FILE" 2>/dev/null)
    if [ $? -ne 0 ]; then
        echo "错误: 无法从 $INFO_FILE 读取群体 $POP 的样本"
        exit 1
    fi
    
    TOTAL_SAMPLES=$(echo "$POP_SAMPLES" | wc -l | awk '{print $1}')
    echo "群体 $POP 的总样本数: $TOTAL_SAMPLES"
    
    if [ $TOTAL_SAMPLES -lt $SAMPLE_PER_POP ]; then
        echo "警告: 样本数不足 ($TOTAL_SAMPLES < $SAMPLE_PER_POP)，使用全部样本"
        SAMPLE_PER_POP=$TOTAL_SAMPLES
    fi
    
    if [ $TOTAL_SAMPLES -eq 0 ]; then
        echo "错误: 群体 $POP 没有样本"
        exit 1
    fi
    
    SAMPLES=$(echo "$POP_SAMPLES" | shuf --random-source="$SEED" 2>/dev/null | head -n $SAMPLE_PER_POP | paste -sd,)
    
    if [ -z "$SAMPLES" ]; then
        echo "错误: 无法获取群体 $POP 的样本"
        exit 1
    fi
    
    echo "抽取样本数: $SAMPLE_PER_POP"
else
    # 合并群体：从两个群体中各抽取一半样本
    SAMPLE_PER_GROUP=$((SAMPLE_PER_POP / 2))
    
    # 检查样本数是否足够
    if [ $SAMPLE_PER_GROUP -lt 1 ]; then
        SAMPLE_PER_GROUP=1
    fi
    
    # 获取第一个群体的样本
    SAMPLES1=$(awk -v pop="$POP1" '$2==pop{print $1}' "$INFO_FILE" 2>/dev/null | shuf --random-source="$SEED" 2>/dev/null | head -n $SAMPLE_PER_GROUP | paste -sd,)
    
    # 获取第二个群体的样本
    SAMPLES2=$(awk -v pop="$POP2" '$2==pop{print $1}' "$INFO_FILE" 2>/dev/null | shuf --random-source="$SEED" 2>/dev/null | head -n $SAMPLE_PER_GROUP | paste -sd,)
    
    # 合并样本
    if [ -n "$SAMPLES1" ] && [ -n "$SAMPLES2" ]; then
        SAMPLES="$SAMPLES1,$SAMPLES2"
    elif [ -n "$SAMPLES1" ]; then
        SAMPLES="$SAMPLES1"
        echo "警告: 只能获取群体 $POP1 的样本"
    elif [ -n "$SAMPLES2" ]; then
        SAMPLES="$SAMPLES2"
        echo "警告: 只能获取群体 $POP2 的样本"
    else
        echo "错误: 无法获取任何样本"
        exit 1
    fi
    
    echo "样本分布: $POP1 ($(echo "$SAMPLES1" | tr ',' '\n' | wc -l)), $POP2 ($(echo "$SAMPLES2" | tr ',' '\n' | wc -l))"
fi

echo "总样本数: $(echo "$SAMPLES" | tr ',' '\n' | wc -l)"

# ===== 输出目录 =====
BASE_DIR=$(pwd)
RESULTS_DIR="${RESULTS_DIR:-results}"
OUTBASE="${BASE_DIR}/${RESULTS_DIR}"
mkdir -p "${OUTBASE}/log"

if [ "$IS_COMBINED" = false ]; then
    POP_NAME="$POP"
    OUTDIR="${OUTBASE}/${POP}"
else
    POP_NAME="$COMBINED_POP"
    OUTDIR="${OUTBASE}/${COMBINED_POP//::/_}"  # 将::替换为_以避免路径问题
fi

mkdir -p "$OUTDIR"/{plots,models,logs,python_scripts}

echo "输出目录: $OUTDIR"

# ===== 检查VCF文件 =====
if [ ! -f "$VCF" ]; then
    echo "错误: VCF文件不存在: $VCF"
    exit 1
fi

# ===== 生成Python脚本 =====
PYTHON_SCRIPT="${OUTDIR}/python_scripts/phlash_${POP_NAME//::/_}.py"

cat > "$PYTHON_SCRIPT" << EOF
#!/usr/bin/env python3
"""
PHLASH分析脚本
分析类型: $ANALYSIS_TYPE
群体: $POP_NAME
样本数: $(echo "$SAMPLES" | tr ',' '\n' | wc -l)
"""

import phlash
import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import traceback
import random
import pickle

def main():
    # 参数设置
    vcf_path = "$VCF"
    mutation_rate = float("$MUTATION_RATE")
    samples = "$SAMPLES".split(',')
    
    # 确定群体名称
    if "$IS_COMBINED" == "true":
        pop = "$COMBINED_POP"
        pop_parts = pop.split("::")
        print(f"分析合并群体: {pop}")
        print(f"包含群体: {pop_parts}")
    else:
        pop = "$POP"
        print(f"分析单个群体: {pop}")
    
    print(f"样本数: {len(samples)}")
    print(f"VCF文件: {vcf_path}")
    print(f"突变率: {mutation_rate}")
    print(f"输出目录: $OUTDIR")
    
    # 染色体长度字典 (GRCh38)
    CHR_LENGTHS = {
        "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
        "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
        "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
        "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
        "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
        "chr21": 46709983, "chr22": 50818468
    }
    
    # 根据MAX_CHROMOSOMES选择染色体
    if int("$MAX_CHROMOSOMES") < 22:
        chromosomes = [f"chr{i}" for i in range(1, int("$MAX_CHROMOSOMES") + 1)]
        print(f"分析前 {len(chromosomes)} 条染色体")
    else:
        chromosomes = [f"chr{i}" for i in range(1, 23)]
    
    print(f"开始PHLASH分析: {datetime.now()}")
    
    # 创建输出目录
    os.makedirs("$OUTDIR/plots", exist_ok=True)
    os.makedirs("$OUTDIR/models", exist_ok=True)
    
    # 参数设置
    WINDOW_SIZE = 1000000  # 1Mb窗口
    MIN_VARIANTS = 10      # 最小变异数阈值
    MAX_CONTIGS_PER_CHROM = 30  # 每条染色体最多使用几个contig
    
    # 构建contigs
    contigs = []
    total_attempts = 0
    successful_attempts = 0
    
    for chrom in chromosomes:
        if chrom not in CHR_LENGTHS:
            print(f"跳过未知染色体: {chrom}")
            continue
        
        chrom_len = CHR_LENGTHS[chrom]
        chrom_contigs = 0
        
        print(f"\n处理染色体: {chrom} (长度: {chrom_len:,} bp)")
        
        # 按窗口遍历染色体
        start_pos = 1
        window_count = 0
        
        while start_pos < chrom_len and chrom_contigs < MAX_CONTIGS_PER_CHROM:
            end_pos = min(start_pos + WINDOW_SIZE - 1, chrom_len)
            region = f"{chrom}:{start_pos}-{end_pos}"
            
            total_attempts += 1
            
            try:
                # 构建contig
                c = phlash.contig(vcf_path, samples=samples, region=region)
                
                # 检查contig是否有有效数据
                if hasattr(c, 'L') and c.L > 0:
                    # 检查变异数量
                    if c.L >= MIN_VARIANTS:
                        contigs.append(c)
                        chrom_contigs += 1
                        successful_attempts += 1
                        print(f"  窗口 {window_count}: 成功构建contig, 变异数: {c.L}")
                    else:
                        print(f"  窗口 {window_count}: 变异数不足 ({c.L} < {MIN_VARIANTS})")
                else:
                    print(f"  窗口 {window_count}: 无有效数据")
                    
            except Exception as e:
                error_msg = str(e)
                if "no samples" in error_msg.lower():
                    print(f"  窗口 {window_count}: 无样本数据")
                elif "no variant" in error_msg.lower():
                    print(f"  窗口 {window_count}: 无变异")
                else:
                    print(f"  窗口 {window_count}: 创建contig失败 - {error_msg[:80]}")
            
            # 移动到下一个窗口
            start_pos += WINDOW_SIZE
            window_count += 1
        
        print(f"染色体 {chrom} 成功构建 {chrom_contigs} 个contig")
    
    print(f"\n窗口尝试总数: {total_attempts}")
    print(f"成功构建contig数: {successful_attempts}")
    print(f"成功率: {successful_attempts/total_attempts*100:.1f}%")
    
    if not contigs:
        raise RuntimeError("没有成功加载任何contig，请检查VCF和样本列表")
    
    print(f"\n总共成功加载 {len(contigs)} 个contig")
    total_variants = sum(c.L for c in contigs)
    print(f"总变异数: {total_variants:,}")
    
    # 检查是否达到最小contig数要求
    if len(contigs) < int("$MIN_CONTIGS"):
        print(f"警告: contig数 ({len(contigs)}) 小于最小要求 ({'$MIN_CONTIGS'})")
        print("可能影响分析结果的可靠性")
    
    # 划分训练和测试数据
    random.shuffle(contigs)
    
    if len(contigs) > 1:
        test_data = contigs[0]
        train_data = contigs[1:]
        print(f"训练数据: {len(train_data)} 个contigs")
        print(f"测试数据: 1 个contig")
    else:
        print("只有一个contig，使用全部数据进行分析")
        train_data = contigs
        test_data = None
    
    # 模型拟合
    print("\n开始拟合模型...")
    try:
        if test_data:
            results = phlash.fit(
                data=train_data,
                test_data=test_data,
                mutation_rate=mutation_rate,
                knots=int("$KNOTS"),
                regularization_penalty=float("$REG_PENALTY"),
                spline_type="$SPLINE_TYPE"
            )
        else:
            results = phlash.fit(
                data=train_data,
                mutation_rate=mutation_rate,
                knots=int("$KNOTS"),
                regularization_penalty=float("$REG_PENALTY"),
                spline_type="$SPLINE_TYPE"
            )
        
        print("模型拟合完成!")
        
        # 保存结果
        with open("$OUTDIR/models/phlash_results.pkl", "wb") as f:
            pickle.dump(results, f)
            
    except Exception as e:
        print(f"模型拟合过程中出错: {e}")
        traceback.print_exc()
        return 1
    
    # 绘图
    times = np.array([dm.eta.t[1:] for dm in results])
    T = np.geomspace(times.min(), times.max(), 1000)
    Nes = np.array([dm.eta(T, Ne=True) for dm in results])
    
    # 绘图1: 所有后验样本
    plt.figure(figsize=(10, 6))
    for ne in Nes:
        plt.plot(T, ne, alpha=0.1, color='blue')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("Time (generations)")
    plt.ylabel("Effective Population Size")
    plt.title(f"{pop} - All Posterior Samples")
    plt.grid(True, alpha=0.3)
    plt.savefig("$OUTDIR/plots/phlash_all_samples.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 绘图2: 中位数和置信区间
    plt.figure(figsize=(10, 6))
    plt.fill_between(
        T,
        np.percentile(Nes, 5, axis=0),
        np.percentile(Nes, 95, axis=0),
        color='lightblue', alpha=0.4, label='90% credible interval'
    )
    plt.plot(T, np.median(Nes, axis=0), color='navy', lw=2, label='Median Ne')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("Time (generations)")
    plt.ylabel("Effective Population Size")
    plt.legend()
    plt.title(f"{pop}: Posterior Median and 90% CI")
    plt.grid(True, alpha=0.3)
    plt.savefig("$OUTDIR/plots/phlash_CI.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 绘图3: 不同百分位数
    plt.figure(figsize=(10, 6))
    percentiles = [10, 25, 50, 75, 90]
    colors = ['lightblue', 'skyblue', 'navy', 'skyblue', 'lightblue']
    
    for i, p in enumerate(percentiles):
        plt.plot(T, np.percentile(Nes, p, axis=0), 
                color=colors[i], 
                alpha=0.6 if i != 2 else 1.0,
                lw=1 if i != 2 else 2,
                label=f'{p}th percentile' if i == 0 or i == len(percentiles)-1 or i == 2 else None)
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("Time (generations)")
    plt.ylabel("Effective Population Size")
    plt.legend()
    plt.title(f"{pop}: Different Percentiles")
    plt.grid(True, alpha=0.3)
    plt.savefig("$OUTDIR/plots/phlash_percentiles.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"分析完成于: {datetime.now()}")
    
    # 生成摘要统计
    print("\n=== 摘要统计 ===")
    print(f"群体: {pop}")
    print(f"分析样本数: {len(samples)}")
    print(f"使用的contig数: {len(contigs)}")
    print(f"总变异数: {total_variants:,}")
    print(f"后验样本数: {len(results)}")
    
    # 计算不同时间点的Ne统计
    time_points = [100, 500, 1000, 5000, 10000, 50000, 100000]
    
    print("\n不同时间点的有效种群大小统计:")
    print(f"{'时间(代)':>10} {'Ne中位数':>12} {'Ne(5%)':>10} {'Ne(95%)':>10}")
    print("-" * 50)
    
    summary_data = []
    for t in time_points:
        # 找到最接近的时间点
        idx = np.argmin(np.abs(T - t))
        ne_median = np.median(Nes[:, idx])
        ne_5th = np.percentile(Nes[:, idx], 5)
        ne_95th = np.percentile(Nes[:, idx], 95)
        
        print(f"{t:10d} {ne_median:12.0f} {ne_5th:10.0f} {ne_95th:10.0f}")
        summary_data.append([t, ne_median, ne_5th, ne_95th])
    
    # 保存摘要统计到文件
    with open("$OUTDIR/models/summary_stats.txt", "w") as f:
        f.write(f"PHLASH Analysis Summary\n")
        f.write(f"=======================\n")
        f.write(f"Population: {pop}\n")
        f.write(f"Analysis type: $ANALYSIS_TYPE\n")
        f.write(f"Sample size: {len(samples)}\n")
        f.write(f"Number of contigs: {len(contigs)}\n")
        f.write(f"Total variants: {total_variants}\n")
        f.write(f"Mutation rate: {mutation_rate}\n")
        f.write(f"Knots: $KNOTS\n")
        f.write(f"Regularization penalty: $REG_PENALTY\n")
        f.write(f"Spline type: $SPLINE_TYPE\n")
        f.write(f"\nEffective Population Size at Different Times:\n")
        f.write(f"{'Time(gens)':>12} {'Ne_median':>12} {'Ne_5th':>12} {'Ne_95th':>12}\n")
        for t, median, p5, p95 in summary_data:
            f.write(f"{t:12d} {median:12.0f} {p5:12.0f} {p95:12.0f}\n")
        f.write(f"\nAnalysis completed: {datetime.now()}\n")
    
    print(f"\n详细结果已保存至: $OUTDIR")
    
    return 0

if __name__ == "__main__":
    exit(main())
EOF

# 设置Python脚本可执行权限
chmod +x "$PYTHON_SCRIPT"

# ===== 运行分析 =====
echo "运行PHLASH分析..."
echo "分析类型: $ANALYSIS_TYPE"
echo "群体: $POP_NAME"
echo "样本数: $(echo "$SAMPLES" | tr ',' '\n' | wc -l)"
echo "Python脚本: $PYTHON_SCRIPT"
echo "输出目录: $OUTDIR"

start_time=$(date +%s)

python "$PYTHON_SCRIPT" > "${OUTDIR}/logs/phlash_run.log" 2>&1

EXIT_STATUS=$?
end_time=$(date +%s)
duration=$((end_time - start_time))

if [ $EXIT_STATUS -eq 0 ]; then
    echo "✓ $POP_NAME 分析成功完成"
    echo "  用时: ${duration}秒"
    echo "  结果保存在: $OUTDIR"
    
    # 检查输出文件
    if [ -f "${OUTDIR}/models/phlash_results.pkl" ]; then
        echo "  模型文件: ✓ 已生成"
    else
        echo "  警告: 模型文件未生成"
    fi
    
    if [ -f "${OUTDIR}/plots/phlash_CI.png" ]; then
        echo "  图表文件: ✓ 已生成"
    else
        echo "  警告: 图表文件未生成"
    fi
else
    echo "✗ $POP_NAME 分析失败 (退出码: $EXIT_STATUS)"
    echo "  用时: ${duration}秒"
    echo "  请检查日志文件: ${OUTDIR}/logs/phlash_run.log"
    
    # 显示日志文件的最后几行
    if [ -f "${OUTDIR}/logs/phlash_run.log" ]; then
        echo "\n=== 日志文件最后20行 ==="
        tail -20 "${OUTDIR}/logs/phlash_run.log"
    fi
    
    exit $EXIT_STATUS
fi

echo "=== PHLASH 分析结束 ==="
echo "完成时间: $(date)"
