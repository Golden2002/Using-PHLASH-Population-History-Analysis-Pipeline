#!/bin/bash
#SBATCH --job-name=PHLASH_SexChrom
#SBATCH --output=/home/litianxing/100My_Jino/114.PHLASH/results/log/PHLASH_sex_%A_%a.log
#SBATCH --error=/home/litianxing/100My_Jino/114.PHLASH/results/log/PHLASH_sex_%A_%a.err.log
#SBATCH --ntasks=1
#SBATCH --partition=batch
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --array=0-2  # 对应性染色体类型：0=X, 1=Y, 2=MT

# ===== 环境设置 =====
source /home/litianxing/100My_Jino/114.PHLASH/software/bin/activate

# ===== 用户可修改参数 =====
POPULATIONS=("Jino" "Han_N" "Tibetan" "Pumi" "Mosuo" "Hani" "Lahu" "Qiang")

# 性染色体VCF路径
VCF_X="/home/litianxing/100My_Jino/101DataPanel/101.1NGS/NGS1926/NGS1926_GRCh38/VCF/200.splite.chr/NGS720_plusRef_1926.VQSR.variants.PASS.SNPs.biallelic.chrX.vcf.gz"
VCF_Y="/share/home/litianxing/100My_Jino/101DataPanel/filtered_NGS1926_chrY.vcf.gz"
VCF_MT="/share/home/litianxing/100My_Jino/101DataPanel/filtered_NGS1926_chrM.vcf.gz"

#MASK_FILE="/share/home/litianxing/100My_Jino/110.Relate/mask/StrictMask/20160622.allChr.mask.bed.gz"
SEED="/home/litianxing/100My_Jino/114.PHLASH/shuf_random_seed.txt"
MUTATION_RATE="1.25e-8"
KNOTS=20
REG_PENALTY=6.0
SPLINE_TYPE="piecewise"

# ===== 男性样本文件 =====
MALE_SAMPLES_FILE="/home/litianxing/100My_Jino/115.Patrilineality/male_samples.txt"

# ===== 当前群体和染色体类型 =====
POPULATIONS=("Jino" "Han_N" "Tibetan" "Pumi" "Mosuo" "Hani" "Lahu" "Qiang")
POP=${POPULATIONS[$((SLURM_ARRAY_TASK_ID / 3))]}  # 8个群体 * 3种性染色体 = 24个任务
CHROM_TYPES=("X" "Y" "MT")
CHROM_TYPE=${CHROM_TYPES[$((SLURM_ARRAY_TASK_ID % 3))]}

echo "分析群体: $POP, 染色体类型: $CHROM_TYPE"

# 提取样本
INFO_FILE="/home/litianxing/00My_Jinuo/01datacollection/infos/filtered_modified_info2.txt"

# 读取男性样本列表
MALE_SAMPLES=$(cat "$MALE_SAMPLES_FILE")

if [ "$CHROM_TYPE" = "Y" ]; then
    # 对于Y染色体，只使用男性样本
    # 获取该群体的所有样本
    POP_SAMPLES=$(awk -v pop="$POP" '$2==pop{print $1}' $INFO_FILE)

    # 取交集：群体样本 ∩ 男性样本
    Y_SAMPLES=$(comm -12 <(echo "$POP_SAMPLES" | sort) <(echo "$MALE_SAMPLES" | sort))

    # 检查是否有足够的男性样本
    Y_SAMPLE_COUNT=$(echo "$Y_SAMPLES" | wc -l)
    echo "群体 $POP 的男性样本数量: $Y_SAMPLE_COUNT"

    if [ $Y_SAMPLE_COUNT -eq 0 ]; then
        echo "错误: 群体 $POP 没有男性样本，无法进行Y染色体分析"
        exit 1
    elif [ $Y_SAMPLE_COUNT -lt 15 ]; then
        echo "警告: 群体 $POP 只有 $Y_SAMPLE_COUNT 个男性样本，将使用全部样本"
        SAMPLES=$(echo "$Y_SAMPLES" | shuf --random-source=$SEED | paste -sd,)
    else
        # 随机选择15个男性样本
        SAMPLES=$(echo "$Y_SAMPLES" | shuf --random-source=$SEED | head -n 15 | paste -sd,)
    fi
else
    # 对于X和MT染色体，使用所有样本
    SAMPLES=$(awk -v pop="$POP" '$2==pop{print $1}' $INFO_FILE | shuf --random-source=$SEED | head -n 15 | paste -sd,)
fi

# ===== 输出目录 =====
OUTBASE="/home/litianxing/100My_Jino/114.PHLASH"
OUTDIR="/home/litianxing/100My_Jino/114.PHLASH/results/${POP}_${CHROM_TYPE}"
mkdir -p $OUTDIR/{plots,models,logs,python_scripts}

# ===== 根据染色体类型设置VCF路径 =====
case $CHROM_TYPE in
    "X") VCF_PATH=$VCF_X ;;
    "Y") VCF_PATH=$VCF_Y ;;
    "MT") VCF_PATH=$VCF_MT ;;
esac

# ===== Python 脚本 =====
PYTHON_SCRIPT="${OUTDIR}/python_scripts/phlash_${POP}_${CHROM_TYPE}.py"
cat > $PYTHON_SCRIPT << EOF
#!/usr/bin/env python3
import phlash, os, numpy as np, matplotlib.pyplot as plt
from datetime import datetime
import traceback, gzip

def main():
    # 根据染色体类型设置参数
    chrom_type = "$CHROM_TYPE"
    pop = "$POP"

    # 设置VCF路径
    if chrom_type == "X":
        vcf_path = "/home/litianxing/100My_Jino/101DataPanel/101.1NGS/NGS1926/NGS1926_GRCh38/VCF/200.splite.chr/NGS720_plusRef_1926.VQSR.variants.PASS.SNPs.biallelic.chrX.vcf.gz"
    elif chrom_type == "Y":
        vcf_path = "/share/home/litianxing/100My_Jino/101DataPanel/filtered_NGS1926_chrY.vcf.gz"
    elif chrom_type == "MT":
        vcf_path = "/share/home/litianxing/100My_Jino/101DataPanel/filtered_NGS1926_chrM.vcf.gz"

    mutation_rate = 1.25e-8
    samples = "$SAMPLES".split(',')

    # ===== 性染色体长度字典 =====
    SEX_CHR_LENGTHS = {
        "chrX": 156040895,  # GRCh38 chrX长度
        "chrY": 57227415,   # GRCh38 chrY长度
        "chrM": 16569       # 线粒体DNA长度
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

    print(f"开始PHLASH分析: {datetime.now()}")
    print(f"群体: {pop}, 染色体类型: {chrom_type}, 样本数: {len(samples)}")
    print(f"前5个样本: {samples[:5]}")

    os.makedirs("${OUTDIR}/plots", exist_ok=True)
    os.makedirs("${OUTDIR}/models", exist_ok=True)

    # ===== 根据染色体类型调整参数 =====
    if chrom_type == "X":
        # X染色体较大，使用较小的窗口
        WINDOW_SIZE = 500000    # 500kb窗口
        MIN_VARIANTS = 5
        MAX_CONTIGS_PER_CHROM = 10
    elif chrom_type == "Y":
        # Y染色体较小，变异较少
        WINDOW_SIZE = 200000    # 200kb窗口
        MIN_VARIANTS = 3
        MAX_CONTIGS_PER_CHROM = 8
    elif chrom_type == "MT":
        # 线粒体很小，不使用窗口
        WINDOW_SIZE = None
        MIN_VARIANTS = 1
        MAX_CONTIGS_PER_CHROM = 1

    # ===== 遍历染色体，划分区域并构建contig =====
    contigs = []

    for chrom in chromosomes:
        if chrom not in SEX_CHR_LENGTHS:
            print(f"跳过未知染色体: {chrom}")
            continue

        chrom_len = SEX_CHR_LENGTHS[chrom]
        chrom_contigs = 0

        print(f"\n处理染色体: {chrom} (长度: {chrom_len})")

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
                        print(f"  成功构建contig, 长度(L): {c.L}")
                    else:
                        print(f"  跳过区域（变异数太少）: {c.L}")
                else:
                    print(f"  跳过区域（无长度信息）")

            except Exception as e:
                print(f"  区域 {region} 创建contig失败: {str(e)[:100]}...")

        else:
            # 对于X和Y染色体，使用窗口化处理
            start_pos = 1
            while start_pos < chrom_len and chrom_contigs < MAX_CONTIGS_PER_CHROM:
                end_pos = min(start_pos + WINDOW_SIZE - 1, chrom_len)
                region = f"{chrom}:{start_pos}-{end_pos}"

                try:
                    print(f"  尝试区域: {region}")
                    c = phlash.contig(vcf_path, samples=samples, region=region)

                    # 检查contig是否有有效数据
                    if hasattr(c, 'L') and c.L > 0:
                        if c.L >= MIN_VARIANTS:
                            contigs.append(c)
                            chrom_contigs += 1
                            print(f"  成功构建contig, 长度(L): {c.L}")
                        else:
                            print(f"  跳过区域（变异数太少）: {c.L}")
                    else:
                        print(f"  跳过区域（无长度信息）")

                except Exception as e:
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
                knots=20,
                regularization_penalty=6.0,
                spline_type="piecewise"
            )
        else:
            results = phlash.fit(
                data=train_data,
                mutation_rate=mutation_rate,
                knots=20,
                regularization_penalty=6.0,
                spline_type="piecewise"
            )

        print("模型拟合完成!")

        # 保存结果
        import pickle
        with open("${OUTDIR}/models/phlash_results.pkl", "wb") as f:
            pickle.dump(results, f)

        # ===== 绘图 =====
        times = np.array([dm.eta.t[1:] for dm in results])
        T = np.geomspace(times.min(), times.max(), 1000)
        Nes = np.array([dm.eta(T, Ne=True) for dm in results])

        plt.figure(figsize=(10,6))
        for ne in Nes:
            plt.plot(T, ne, alpha=0.1, color='blue')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel("time (in generations)")
        plt.ylabel("Effective Population Size")
        plt.title(f"{pop} {chrom_type} - all posterior samples")
        plt.savefig("${OUTDIR}/plots/phlash_all_samples.png", dpi=300, bbox_inches='tight')

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
        plt.title(f"{pop} {chrom_type}: Posterior Median and 90% CI")
        plt.savefig("${OUTDIR}/plots/phlash_CI.png", dpi=300, bbox_inches='tight')

        print(f"分析完成于: {datetime.now()}")

    except Exception as e:
        print(f"模型拟合过程中出错: {e}")
        traceback.print_exc()
        return 1

    return 0

if __name__ == "__main__":
    exit(main())
EOF

# ===== 运行 =====
echo "运行PHLASH分析: 群体 $POP, 染色体 $CHROM_TYPE"
python $PYTHON_SCRIPT > ${OUTDIR}/logs/phlash_run.log 2>&1

if [ $? -eq 0 ]; then
    echo "群体 $POP 染色体 $CHROM_TYPE 分析完成"
else
    echo "群体 $POP 染色体 $CHROM_TYPE 分析失败，请检查日志: ${OUTDIR}/logs/phlash_run.log"
    exit 1
fi