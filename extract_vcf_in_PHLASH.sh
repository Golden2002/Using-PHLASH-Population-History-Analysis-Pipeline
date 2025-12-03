#!/bin/bash
#SBATCH --job-name=VCF_Extraction
#SBATCH --output=logs/VCF_Extraction.log
#SBATCH --error=logs/VCF_Extraction.err.log
#SBATCH --ntasks=1
#SBATCH --partition=batch
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1

# ===== 环境设置 =====
# 请根据实际环境修改虚拟环境路径
source /path/to/your/conda_or_venv/bin/activate

# ===== 用户可修改参数 =====
POPULATIONS=("POP1" "POP2" "POP3")     # 目标群体列表，可增加或修改
VCF="/path/to/input/your_data.vcf.gz" # 输入VCF文件
MASK_FILE="/path/to/mask_file.bed.gz" # mask区域文件，bgzip压缩
INFO_FILE="/path/to/sample_info.txt"  # 样本信息文件，包含: 样本ID 群体名

# ===== 输出目录 =====
OUTBASE="./output"
mkdir -p $OUTBASE
mkdir -p logs

ALL_SAMPLES="$OUTBASE/ALL_TARGET_SAMPLES.txt"
VCF_STEP1_TMP="$OUTBASE/step1.tmp.vcf.gz"
VCF_STEP1="$OUTBASE/step1.filtered.vcf.gz"
VCF_FINAL="$OUTBASE/final_filtered.vcf.gz"

# ===== Step 1: 收集群体样本 =====
echo ">>> Step 1: 收集所有群体样本"
rm -f $ALL_SAMPLES
for p in "${POPULATIONS[@]}"; do
    awk -v pop="$p" '$2==pop{print $1}' $INFO_FILE >> $ALL_SAMPLES
done
echo "总样本数：$(wc -l < $ALL_SAMPLES)"

# ===== Step 2a: 样本过滤 =====
echo ">>> Step 2a: 过滤样本"
bcftools view \
    -S $ALL_SAMPLES \
    -m2 -M2 -v snps \
    $VCF \
    -Oz -o $VCF_STEP1_TMP
bcftools index -t $VCF_STEP1_TMP

# ===== Step 2b: MAC过滤 =====
echo ">>> Step 2b: 过滤次要等位基因数(MAC)"
vcftools \
    --gzvcf $VCF_STEP1_TMP \
    --mac 2 \
    --recode \
    --recode-INFO-all \
    --stdout \
    | bgzip -c > $VCF_STEP1
tabix -p vcf $VCF_STEP1

# 清理临时文件
rm -f $VCF_STEP1_TMP ${VCF_STEP1_TMP}.tbi

# ===== Step 3: 去除 mask 区域 =====
echo ">>> Step 3: 去除 mask 区域"
bcftools view \
    -T ^$MASK_FILE $VCF_STEP1 \
    -Oz -o $VCF_FINAL
bcftools index -t $VCF_FINAL

echo ">>> VCF预处理完成：$VCF_FINAL"
