#!/bin/bash
#SBATCH --job-name=VCF_Extraction
#SBATCH --output=/home/litianxing/100My_Jino/114.PHLASH/results/logs/VCF_Extraction.log
#SBATCH --error=/home/litianxing/100My_Jino/114.PHLASH/results/logs/VCF_Extraction.err.log
#SBATCH --ntasks=1
#SBATCH --partition=batch
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1

# ===== 环境设置 =====
source /home/litianxing/100My_Jino/114.PHLASH/software/bin/activate

# ===== 用户可修改参数 =====
POPULATIONS=("Jino" "Han_N" "Tibetan" "Pumi" "Mosuo" "Hani" "Lahu" "Qiang")
VCF="/share/home/litianxing/100My_Jino/107.IBD/data/NGS.phased.vcf.gz"
MASK_FILE="/share/home/litianxing/100My_Jino/110.Relate/mask/StrictMask/20160622.allChr.mask.bed.gz"
INFO_FILE="/home/litianxing/00My_Jinuo/01datacollection/infos/filtered_modified_info2.txt"

# ===== 输出目录 =====
OUTBASE="/home/litianxing/100My_Jino/114.PHLASH/data"
ALL_SAMPLES="$OUTBASE/ALL_TARGET_SAMPLES.txt"
VCF_STEP1="$OUTBASE/clean.step1.samples_maf.vcf.gz"
VCF_FINAL="$OUTBASE/clean.ALL_TARGETS.vcf.gz"
VCF_STEP1_TMP="$OUTBASE/clean.step1.tmp.vcf.gz"

echo ">>> Step 1: 收集所有群体样本"
rm -f $ALL_SAMPLES
for p in "${POPULATIONS[@]}"; do
    awk -v pop="$p" '$2==pop{print $1}' $INFO_FILE >> $ALL_SAMPLES
done

echo "总样本数：$(wc -l < $ALL_SAMPLES)"

echo ">>> Step 2a: 过滤样本"
bcftools view \
    -S $ALL_SAMPLES \
    -m2 -M2 -v snps \
    $VCF \
    -Oz -o $VCF_STEP1_TMP
bcftools index -t $VCF_STEP1_TMP

echo ">>> Step 2b: 过滤MAC"
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

echo ">>> Step 3: 去除 mask 区域"
bcftools view \
    -T ^$MASK_FILE $VCF_STEP1 \
    -Oz -o $VCF_FINAL
bcftools index -t $VCF_FINAL

echo ">>> VCF预处理完成：$VCF_FINAL"