#!/bin/bash
#SBATCH --job-name=PHLASH_SexChrom
#SBATCH --output=results/log/PHLASH_sex_%A_%a.log
#SBATCH --error=results/log/PHLASH_sex_%A_%a.err.log
#SBATCH --ntasks=1
#SBATCH --partition=batch
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --array=0-2  # 对应性染色体类型：0=X, 1=Y, 2=MT

# ===== 配置加载 =====
# 从外部配置文件加载敏感信息和路径
if [ -f "config.sh" ]; then
    source config.sh
else
    echo "错误: 找不到配置文件 config.sh"
    echo "请创建 config.sh 并设置以下变量:"
    echo "  VCF_X, VCF_Y, VCF_MT, INFO_FILE, MALE_SAMPLES_FILE"
    exit 1
fi

# ===== 环境设置 =====
# 激活conda环境或设置软件路径
if [ -n "$PHLASH_ENV" ]; then
    source "$PHLASH_ENV"
fi

# ===== 用户可修改参数 =====
# 群体列表 - 根据实际研究修改
POPULATIONS=("Pop1" "Pop2" "Pop3" "Pop4" "Pop5" "Pop6" "Pop7" "Pop8")

# 检查必要的变量是否设置
if [ -z "$VCF_X" ] || [ -z "$VCF_Y" ] || [ -z "$VCF_MT" ]; then
    echo "错误: VCF路径未在config.sh中设置"
    exit 1
fi

# 默认参数
MUTATION_RATE=${MUTATION_RATE:-"1.25e-8"}
KNOTS=${KNOTS:-20}
REG_PENALTY=${REG_PENALTY:-6.0}
SPLINE_TYPE=${SPLINE_TYPE:-"piecewise"}
SEED=${SEED:-"shuf_random_seed.txt"}
SAMPLE_SIZE=${SAMPLE_SIZE:-15}  # 每个群体使用的样本数
MIN_MALE_SAMPLES=${MIN_MALE_SAMPLES:-15}  # Y染色体分析最小男性样本数

# ===== 输出目录结构 =====
BASE_DIR=$(pwd)
RESULTS_DIR="${RESULTS_DIR:-results}"
OUTBASE="${BASE_DIR}/${RESULTS_DIR}"
mkdir -p "${OUTBASE}/log"

# ===== 当前群体和染色体类型 =====
# 计算群体和染色体类型
POP_COUNT=${#POPULATIONS[@]}
POP_INDEX=$((SLURM_ARRAY_TASK_ID / 3))
CHROM_INDEX=$((SLURM_ARRAY_TASK_ID % 3))

# 检查索引是否有效
if [ $POP_INDEX -ge $POP_COUNT ]; then
    echo "警告: 任务ID $SLURM_ARRAY_TASK_ID 超出范围"
    echo "      最大任务数应为 $((POP_COUNT * 3 - 1))"
    exit 0  # 正常退出，不报错
fi

POP=${POPULATIONS[$POP_INDEX]}
CHROM_TYPES=("X" "Y" "MT")
CHROM_TYPE=${CHROM_TYPES[$CHROM_INDEX]}

echo "=== PHLASH 分析开始 ==="
echo "时间: $(date)"
echo "任务ID: $SLURM_ARRAY_TASK_ID"
echo "分析群体: $POP (索引: $POP_INDEX)"
echo "染色体类型: $CHROM_TYPE (索引: $CHROM_INDEX)"
echo "工作目录: $BASE_DIR"

# ===== 样本提取 =====
if [ "$CHROM_TYPE" = "Y" ]; then
    # 对于Y染色体，只使用男性样本
    if [ ! -f "$MALE_SAMPLES_FILE" ]; then
        echo "错误: 男性样本文件不存在: $MALE_SAMPLES_FILE"
        exit 1
    fi
    
    # 获取该群体的所有样本
    if [ ! -f "$INFO_FILE" ]; then
        echo "错误: 样本信息文件不存在: $INFO_FILE"
        exit 1
    fi
    
    POP_SAMPLES=$(awk -v pop="$POP" '$2==pop{print $1}' "$INFO_FILE" 2>/dev/null)
    if [ $? -ne 0 ]; then
        echo "错误: 无法从 $INFO_FILE 读取群体 $POP 的样本"
        exit 1
    fi
    
    # 读取男性样本
    MALE_SAMPLES=$(cat "$MALE_SAMPLES_FILE")
    
    # 取交集：群体样本 ∩ 男性样本
    Y_SAMPLES=$(comm -12 <(echo "$POP_SAMPLES" | sort) <(echo "$MALE_SAMPLES" | sort))
    Y_SAMPLE_COUNT=$(echo "$Y_SAMPLES" | wc -l | awk '{print $1}')
    
    echo "群体 $POP 的男性样本数量: $Y_SAMPLE_COUNT"
    
    if [ $Y_SAMPLE_COUNT -eq 0 ]; then
        echo "错误: 群体 $POP 没有男性样本，无法进行Y染色体分析"
        exit 1
    elif [ $Y_SAMPLE_COUNT -lt $MIN_MALE_SAMPLES ]; then
        echo "警告: 群体 $POP 只有 $Y_SAMPLE_COUNT 个男性样本，将使用全部样本"
        SAMPLES=$(echo "$Y_SAMPLES" | shuf --random-source="$SEED" | paste -sd,)
    else
        SAMPLES=$(echo "$Y_SAMPLES" | shuf --random-source="$SEED" | head -n $SAMPLE_SIZE | paste -sd,)
    fi
else
    # 对于X和MT染色体，使用所有样本
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
    echo "群体 $POP 的总样本数量: $TOTAL_SAMPLES"
    
    if [ $TOTAL_SAMPLES -lt $SAMPLE_SIZE ]; then
        echo "警告: 群体 $POP 只有 $TOTAL_SAMPLES 个样本，将使用全部样本"
        SAMPLE_SIZE=$TOTAL_SAMPLES
    fi
    
    SAMPLES=$(echo "$POP_SAMPLES" | shuf --random-source="$SEED" | head -n $SAMPLE_SIZE | paste -sd,)
fi

# 检查是否成功获取样本
if [ -z "$SAMPLES" ]; then
    echo "错误: 无法获取群体 $POP 的样本"
    exit 1
fi

echo "使用的样本数: $(echo "$SAMPLES" | tr ',' '\n' | wc -l)"

# ===== 创建输出目录 =====
OUTDIR="${OUTBASE}/${POP}_${CHROM_TYPE}"
mkdir -p "$OUTDIR"/{plots,models,logs,python_scripts}

# ===== 根据染色体类型设置VCF路径 =====
case $CHROM_TYPE in
    "X") VCF_PATH="$VCF_X" ;;
    "Y") VCF_PATH="$VCF_Y" ;;
    "MT") VCF_PATH="$VCF_MT" ;;
esac

# 检查VCF文件是否存在
if [ ! -f "$VCF_PATH" ]; then
    echo "错误: VCF文件不存在: $VCF_PATH"
    exit 1
fi

# ===== 生成Python脚本 =====
PYTHON_SCRIPT="${OUTDIR}/python_scripts/phlash_${POP}_${CHROM_TYPE}.py"

cat > "$PYTHON_SCRIPT" << EOF
#!/usr/bin/env python3
"""
PHLASH分析脚本
群体: $POP
染色体类型: $CHROM_TYPE
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
    chrom_type = "$CHROM_TYPE"
    pop = "$POP"
    vcf_path = "$VCF_PATH"
    mutation_rate = float("$MUTATION_RATE")
    samples = "$SAMPLES".split(',')
    
    # 性染色体长度 (GRCh38)
    SEX_CHR_LENGTHS = {
        "chrX": 156040895,
        "chrY": 57227415,
        "chrM": 16569
    }
    
    # 根据染色体类型选择染色体名称
    if chrom_type == "X":
        chromosomes = ["chrX"]
        chrom_name = "chrX"
    elif chrom_type == "Y":
        chromosomes = ["chrY"]
        chrom_name = "chrY"
    elif chrom_type == "MT":
        chromosomes = ["chrM"]
        chrom_name = "chrM"
    
    print(f"=== PHLASH 分析开始 ===")
    print(f"时间: {datetime.now()}")
    print(f"群体: {pop}")
    print(f"染色体类型: {chrom_type}")
    print(f"样本数: {len(samples)}")
    print(f"VCF文件: {vcf_path}")
    print(f"突变率: {mutation_rate}")
    print(f"输出目录: $OUTDIR")
    
    # 创建输出目录
    os.makedirs("$OUTDIR/plots", exist_ok=True)
    os.makedirs("$OUTDIR/models", exist_ok=True)
    
    # 根据染色体类型调整参数
    if chrom_type == "X":
        WINDOW_SIZE = 500000    # 500kb窗口
        MIN_VARIANTS = 5
        MAX_CONTIGS_PER_CHROM = 10
    elif chrom_type == "Y":
        WINDOW_SIZE = 200000    # 200kb窗口
        MIN_VARIANTS = 3
        MAX_CONTIGS_PER_CHROM = 8
    elif chrom_type == "MT":
        WINDOW_SIZE = None
        MIN_VARIANTS = 1
        MAX_CONTIGS_PER_CHROM = 1
    
    # 构建contigs
    contigs = []
    
    for chrom in chromosomes:
        if chrom not in SEX_CHR_LENGTHS:
            print(f"跳过未知染色体: {chrom}")
            continue
        
        chrom_len = SEX_CHR_LENGTHS[chrom]
        chrom_contigs = 0
        
        print(f"\n处理染色体: {chrom} (长度: {chrom_len:,} bp)")
        
        if chrom_type == "MT":
            # 对于MT，直接使用整个染色体
            region = f"{chrom}:1-{chrom_len}"
            try:
                print(f"  尝试区域: {region}")
                c = phlash.contig(vcf_path, samples=samples, region=region)
                
                if hasattr(c, 'L') and c.L > 0:
                    if c.L >= MIN_VARIANTS:
                        contigs.append(c)
                        chrom_contigs += 1
                        print(f"  成功构建contig, 变异数: {c.L}")
                    else:
                        print(f"  跳过区域（变异数太少）: {c.L}")
                else:
                    print(f"  跳过区域（无长度信息）")
                    
            except Exception as e:
                print(f"  区域 {region} 创建contig失败: {str(e)[:100]}")
        else:
            # 对于X和Y染色体，使用窗口化处理
            start_pos = 1
            while start_pos < chrom_len and chrom_contigs < MAX_CONTIGS_PER_CHROM:
                end_pos = min(start_pos + WINDOW_SIZE - 1, chrom_len)
                region = f"{chrom}:{start_pos}-{end_pos}"
                
                try:
                    print(f"  尝试区域: {region}")
                    c = phlash.contig(vcf_path, samples=samples, region=region)
                    
                    if hasattr(c, 'L') and c.L > 0:
                        if c.L >= MIN_VARIANTS:
                            contigs.append(c)
                            chrom_contigs += 1
                            print(f"  成功构建contig, 变异数: {c.L}")
                        else:
                            print(f"  跳过区域（变异数太少）: {c.L}")
                    else:
                        print(f"  跳过区域（无长度信息）")
                        
                except Exception as e:
                    print(f"  区域 {region} 创建contig失败: {str(e)[:100]}")
                
                # 移动到下一个窗口
                start_pos += WINDOW_SIZE
        
        print(f"染色体 {chrom} 成功构建 {chrom_contigs} 个contig")
    
    if not contigs:
        raise RuntimeError("没有成功加载任何contig，请检查VCF和样本列表")
    
    print(f"\n总共成功加载 {len(contigs)} 个contig")
    print(f"总变异数: {sum(c.L for c in contigs):,}")
    
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
                knots=$KNOTS,
                regularization_penalty=$REG_PENALTY,
                spline_type="$SPLINE_TYPE"
            )
        else:
            results = phlash.fit(
                data=train_data,
                mutation_rate=mutation_rate,
                knots=$KNOTS,
                regularization_penalty=$REG_PENALTY,
                spline_type="$SPLINE_TYPE"
            )
        
        print("模型拟合完成!")
        
        # 保存结果
        with open("$OUTDIR/models/phlash_results.pkl", "wb") as f:
            pickle.dump(results, f)
        
        # 绘图
        times = np.array([dm.eta.t[1:] for dm in results])
        T = np.geomspace(times.min(), times.max(), 1000)
        Nes = np.array([dm.eta(T, Ne=True) for dm in results])
        
        # 所有后验样本图
        plt.figure(figsize=(10, 6))
        for ne in Nes:
            plt.plot(T, ne, alpha=0.1, color='blue')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel("Time (generations)")
        plt.ylabel("Effective Population Size")
        plt.title(f"{pop} {chrom_type} - All Posterior Samples")
        plt.grid(True, alpha=0.3)
        plt.savefig("$OUTDIR/plots/phlash_all_samples.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # 中位数和可信区间图
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
        plt.title(f"{pop} {chrom_type}: Posterior Median and 90% CI")
        plt.grid(True, alpha=0.3)
        plt.savefig("$OUTDIR/plots/phlash_CI.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # 保存总结统计信息
        with open("$OUTDIR/models/summary_stats.txt", "w") as f:
            f.write(f"PHLASH Analysis Summary\n")
            f.write(f"=======================\n")
            f.write(f"Population: {pop}\n")
            f.write(f"Chromosome: {chrom_type}\n")
            f.write(f"Sample size: {len(samples)}\n")
            f.write(f"Number of contigs: {len(contigs)}\n")
            f.write(f"Total variants: {sum(c.L for c in contigs)}\n")
            f.write(f"Mutation rate: {mutation_rate}\n")
            f.write(f"Knots: $KNOTS\n")
            f.write(f"Regularization penalty: $REG_PENALTY\n")
            f.write(f"Spline type: $SPLINE_TYPE\n")
            f.write(f"Analysis completed: {datetime.now()}\n")
        
        print(f"\n分析完成于: {datetime.now()}")
        print(f"结果保存至: $OUTDIR")
        
    except Exception as e:
        print(f"模型拟合过程中出错: {e}")
        traceback.print_exc()
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
EOF

# 设置Python脚本可执行权限
chmod +x "$PYTHON_SCRIPT"

# ===== 运行分析 =====
echo "运行PHLASH分析..."
echo "群体: $POP"
echo "染色体: $CHROM_TYPE"
echo "Python脚本: $PYTHON_SCRIPT"
echo "输出目录: $OUTDIR"

python "$PYTHON_SCRIPT" > "${OUTDIR}/logs/phlash_run.log" 2>&1

EXIT_STATUS=$?

if [ $EXIT_STATUS -eq 0 ]; then
    echo "✓ 群体 $POP 染色体 $CHROM_TYPE 分析成功完成"
    echo "结果保存在: $OUTDIR"
else
    echo "✗ 群体 $POP 染色体 $CHROM_TYPE 分析失败 (退出码: $EXIT_STATUS)"
    echo "请检查日志文件: ${OUTDIR}/logs/phlash_run.log"
    exit $EXIT_STATUS
fi

echo "=== PHLASH 分析结束 ==="
echo "完成时间: $(date)"
