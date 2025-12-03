#!/bin/bash
#SBATCH --job-name=PHLASH_Analysis
#SBATCH --output=/home/litianxing/100My_Jino/114.PHLASH/results/log/PHLASH_%A_%a.log
#SBATCH --error=/home/litianxing/100My_Jino/114.PHLASH/results/log/PHLASH_%A_%a.err.log
#SBATCH --ntasks=1
#SBATCH --partition=batch
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --array=3-3  # 更新索引：0-2为单个群体，3为合并分析

# ===== 环境设置 =====
source /home/litianxing/100My_Jino/114.PHLASH/software/bin/activate

# ===== 用户可修改参数 =====
POPULATIONS=("Jino" "Han_N" "Tibetan")  # 单个群体
COMBINED_POPULATIONS=("Jino::Han_N")    # 合并群体
INFO_FILE="/home/litianxing/00My_Jinuo/01datacollection/infos/filtered_modified_info2.txt"
VCF="/share/home/litianxing/100My_Jino/107.IBD/data/NGS.phased.vcf.gz"
MASK_FILE="/share/home/litianxing/100My_Jino/110.Relate/mask/StrictMask/20160622.allChr.mask.bed.gz"
SEED="/home/litianxing/100My_Jino/114.PHLASH/shuf_random_seed.txt"
MUTATION_RATE="1.25e-8"
KNOTS=20
REG_PENALTY=6.0
SPLINE_TYPE="piecewise"
SAMPLE_PER_POP=10  # 每个群体抽取的样本数

# ===== 根据任务索引决定分析类型 =====
if [ $SLURM_ARRAY_TASK_ID -lt ${#POPULATIONS[@]} ]; then
  # 单个群体分析
  POP=${POPULATIONS[$SLURM_ARRAY_TASK_ID]}
  IS_COMBINED=false
  echo "分析单个群体: $POP"
else
  # 合并群体分析
  COMBINED_IDX=$((SLURM_ARRAY_TASK_ID - ${#POPULATIONS[@]}))
  COMBINED_POP=${COMBINED_POPULATIONS[$COMBINED_IDX]}
  IS_COMBINED=true
  echo "分析合并群体: $COMBINED_POP"

  # 解析合并群体名称（例如"Jino::Han_N"）
  IFS='::' read -ra POP_ARRAY <<< "$COMBINED_POP"
  POP1=${POP_ARRAY[0]}
  POP2=${POP_ARRAY[1]}
  echo "合并群体包含: $POP1 和 $POP2"
fi

# ===== 提取样本 =====
if [ "$IS_COMBINED" = false ]; then
  # 单个群体：抽取指定数量的样本
  SAMPLES=$(awk -v pop="$POP" '$2==pop{print $1}' $INFO_FILE | shuf --random-source=$SEED | head -n $SAMPLE_PER_POP | paste -sd,)
  echo "抽取 $SAMPLE_PER_POP 个样本: $SAMPLES"
else
  # 合并群体：从两个群体中各抽取一半样本
  SAMPLES1=$(awk -v pop="$POP1" '$2==pop{print $1}' $INFO_FILE | shuf --random-source=$SEED | head -n $((SAMPLE_PER_POP/2)) | paste -sd,)
  SAMPLES2=$(awk -v pop="$POP2" '$2==pop{print $1}' $INFO_FILE | shuf --random-source=$SEED | head -n $((SAMPLE_PER_POP/2)) | paste -sd,)

  # 合并两个群体的样本
  if [ -n "$SAMPLES1" ] && [ -n "$SAMPLES2" ]; then
    SAMPLES="$SAMPLES1,$SAMPLES2"
  elif [ -n "$SAMPLES1" ]; then
    SAMPLES="$SAMPLES1"
  elif [ -n "$SAMPLES2" ]; then
    SAMPLES="$SAMPLES2"
  fi

  echo "抽取合并样本: $SAMPLES1 和 $SAMPLES2"
  echo "总样本数: $(echo $SAMPLES | tr ',' '\n' | wc -l)"
fi

# ===== 输出目录 =====
OUTBASE="/home/litianxing/100My_Jino/114.PHLASH"

if [ "$IS_COMBINED" = false ]; then
  OUTDIR="/home/litianxing/100My_Jino/114.PHLASH/results/$POP"
  POP_NAME="$POP"
else
  # 合并群体的输出目录
  OUTDIR="/home/litianxing/100My_Jino/114.PHLASH/results/$COMBINED_POP"
  POP_NAME="$COMBINED_POP"
fi

mkdir -p $OUTDIR/{plots,models,logs,python_scripts}

# ===== Python 脚本 =====
PYTHON_SCRIPT="${OUTDIR}/python_scripts/phlash_${POP_NAME}.py"
cat > $PYTHON_SCRIPT << EOF
#!/usr/bin/env python3
import phlash, os, numpy as np, matplotlib.pyplot as plt
from datetime import datetime
import traceback, gzip

def main():
    vcf_path = "/share/home/litianxing/100My_Jino/107.IBD/data/NGS.phased.vcf.gz"
    mutation_rate = ${MUTATION_RATE}
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
    print(f"样本列表: {samples}")

    # ===== 染色体长度字典 =====
    CHR_LENGTHS = {
        "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
        "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
        "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
        "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
        "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
        "chr21": 46709983, "chr22": 50818468
    }

    print(f"开始PHLASH分析: {datetime.now()}")
    os.makedirs("${OUTDIR}/plots", exist_ok=True)
    os.makedirs("${OUTDIR}/models", exist_ok=True)

    # ===== 参数设置 =====
    WINDOW_SIZE = 1000000  # 1Mb窗口
    MIN_VARIANTS = 10      # 最小变异数阈值
    MAX_CONTIGS_PER_CHROM = 30  # 每条染色体最多使用几个contig

    # ===== 遍历染色体，划分区域并构建contig =====
    contigs = []

    # 选择要分析的染色体（为了节省时间，可以先测试少量染色体）
    chromosomes = [f"chr{i}" for i in range(1, 23)]  # 全部染色体
    # chromosomes = ["chr1", "chr2"]  # 测试用少量染色体

    for chrom in chromosomes:
        if chrom not in CHR_LENGTHS:
            print(f"跳过未知染色体: {chrom}")
            continue

        chrom_len = CHR_LENGTHS[chrom]
        chrom_contigs = 0

        print(f"\n处理染色体: {chrom} (长度: {chrom_len})")

        # 按窗口遍历染色体
        start_pos = 1
        while start_pos < chrom_len and chrom_contigs < MAX_CONTIGS_PER_CHROM:
            end_pos = min(start_pos + WINDOW_SIZE - 1, chrom_len)
            region = f"{chrom}:{start_pos}-{end_pos}"

            try:
                print(f"  尝试区域: {region}")

                # 构建contig
                c = phlash.contig(vcf_path, samples=samples, region=region)

                # 检查contig是否有有效数据
                if hasattr(c, 'L') and c.L > 0:
                    # 检查变异数量（通过L属性判断）
                    if c.L >= MIN_VARIANTS:
                        contigs.append(c)
                        chrom_contigs += 1
                        print(f"  成功构建contig, 长度(L): {c.L}")
                    else:
                        print(f"  跳过区域（变异数太少）: {c.L}")
                else:
                    print(f"  跳过区域（无长度信息）")

            except Exception as e:
                # 如果出现错误，可能是该区域没有数据或格式问题
                print(f"  区域 {region} 创建contig失败: {str(e)[:100]}...")

            # 移动到下一个窗口
            start_pos += WINDOW_SIZE

        print(f"染色体 {chrom} 成功构建 {chrom_contigs} 个contig")

    if not contigs:
        raise RuntimeError("没有成功加载任何 contig，请检查VCF和样本列表")

    print(f"\n总共成功加载 {len(contigs)} 个 contig")
    print(f"总长度: {sum(c.L for c in contigs)} bp")

    # ===== 划分训练和测试数据 =====
    import random
    random.shuffle(contigs)

    if len(contigs) > 1:
        test_data = contigs[0]  # 单个contig作为测试数据
        train_data = contigs[1:]  # 其余作为训练数据
        print(f"训练数据: {len(train_data)} 个contigs")
        print(f"测试数据: 1 个contig")
    else:
        print("只有一个contig，使用全部数据进行分析")
        train_data = contigs
        test_data = None

    # ===== 模型拟合 =====
    print("开始拟合模型...")
    try:
        if test_data:
            results = phlash.fit(
                data=train_data,
                test_data=test_data,
                mutation_rate=mutation_rate,
                knots=${KNOTS},
                regularization_penalty=${REG_PENALTY},
                spline_type="${SPLINE_TYPE}"
            )
        else:
            results = phlash.fit(
                data=train_data,
                mutation_rate=mutation_rate,
                knots=${KNOTS},
                regularization_penalty=${REG_PENALTY},
                spline_type="${SPLINE_TYPE}"
            )

        print("模型拟合完成!")

        # 保存结果
        import pickle
        with open("${OUTDIR}/models/phlash_results.pkl", "wb") as f:
            pickle.dump(results, f)
    except Exception as e:
        print(f"模型拟合过程中出错: {e}")
        traceback.print_exc()
        return 1

    # ===== 绘图 =====
    times = np.array([dm.eta.t[1:] for dm in results])
    T = np.geomspace(times.min(), times.max(), 1000)
    Nes = np.array([dm.eta(T, Ne=True) for dm in results])

    # 绘图1: 所有后验样本
    plt.figure(figsize=(10,6))
    for ne in Nes:
        plt.plot(T, ne, alpha=0.1, color='blue')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("time (in generations)")
    plt.ylabel("Effective Population Size")
    plt.title(f"{pop} - all posterior samples")
    plt.savefig("${OUTDIR}/plots/phlash_all_samples.png", dpi=300, bbox_inches='tight')

    # 绘图2: 中位数和置信区间
    plt.figure(figsize=(10,6))
    plt.fill_between(
        T,
        np.percentile(Nes, 5, axis=0),
        np.percentile(Nes, 95, axis=0),
        color='lightblue', alpha=0.4, label='90% credible interval'
    )
    plt.plot(T, np.median(Nes, axis=0), color='navy', lw=2, label='Median Ne')
    plt.xscale('log'); plt.yscale('log')
    plt.xlabel("Time (generations)")
    plt.ylabel("Effective Population Size")
    plt.legend()
    plt.title(f"{pop}: Posterior Median and 90% CI")
    plt.savefig("${OUTDIR}/plots/phlash_CI.png", dpi=300, bbox_inches='tight')

    print(f"分析完成于: {datetime.now()}")

    # ===== 生成摘要统计 =====
    print("\n=== 摘要统计 ===")
    print(f"群体: {pop}")
    print(f"分析样本数: {len(samples)}")
    print(f"使用的contig数: {len(contigs)}")
    print(f"总序列长度: {sum(c.L for c in contigs):,} bp")
    print(f"后验样本数: {len(results)}")

    # 计算不同时间点的Ne中位数
    time_points = [100, 1000, 10000, 100000]
    ne_median = np.median(Nes, axis=0)

    print("\n不同时间点的有效种群大小中位数:")
    for t in time_points:
        # 找到最接近的时间点
        idx = np.argmin(np.abs(T - t))
        print(f"  {t:6d} 代: Ne ≈ {ne_median[idx]:.0f}")

if __name__ == "__main__":
    exit(main())
EOF

# ===== 运行 =====
echo "运行PHLASH分析..."
python $PYTHON_SCRIPT > ${OUTDIR}/logs/phlash_run.log 2>&1

if [ $? -eq 0 ]; then
    echo "$POP_NAME 分析完成"
else
    echo "$POP_NAME 分析失败，请检查日志: ${OUTDIR}/logs/phlash_run.log"
    exit 1
fi