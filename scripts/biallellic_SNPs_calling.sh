#!/bin/bash
#SBATCH --job-name=variant_calling_biallelicSNPs
#SBATCH -p long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128000
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=v.nascimento@uniandes.edu.co
#SBATCH --mail-type=ALL
#SBATCH -o variant_calling_%j.out
#SBATCH -e variant_calling_%j.err

################### modules
module load bwa-mem2/2.2.1
module load samtools/1.16.1
module load bcftools/1.16
module load htslib/1.16

################### variables
THREADS=32
REF=hifiasm_Sample_X.p_ctg.fasta
SAMPLES=(Sample_10 Sample_7 Sample_27 Sample_21 Sample_26 Alim_37 Agly_50)

################### index reference
bwa-mem2 index $REF
samtools faidx $REF

################### mapping and deduplication
for S in "${SAMPLES[@]}"; do
  bwa-mem2 mem -t $THREADS -R "@RG\tID:${S}\tSM:${S}\tPL:ILLUMINA" \
    $REF ${S}_R1.paired.fastq.gz ${S}_R2.paired.fastq.gz \
    | samtools sort -@ $THREADS -o ${S}.sorted.bam -
  samtools markdup -@ $THREADS ${S}.sorted.bam ${S}.dedup.bam
  samtools index ${S}.dedup.bam
  samtools flagstat ${S}.dedup.bam > ${S}.flagstat.txt
done

################### bam list
BAMLIST=bam_list.txt
printf "%s\n" "${SAMPLES[@]/%/.dedup.bam}" > $BAMLIST

################### joint calling
# I prefer bcftools over GATK here, no known-sites VCF for BQSR and lighter on the HPC
bcftools mpileup -f $REF -b $BAMLIST -a AD,DP --threads $THREADS -Ou \
  | bcftools call -mv --threads $THREADS -Oz -o Sample_X.raw.vcf.gz
tabix -p vcf Sample_X.raw.vcf.gz

################### filtering
# biallelic SNPs only, present in at least 80% of samples, QUAL >= 30
bcftools view -m2 -M2 -v snps Sample_X.raw.vcf.gz -Oz \
  | bcftools view -e 'F_MISSING > 0.2' -Oz \
  | bcftools filter -e 'QUAL<30' -Oz -o Sample_X.filtered.vcf.gz
tabix -p vcf Sample_X.filtered.vcf.gz

################### SNP statistics
# overall counts, ts/tv ratio, quality and depth distributions, one number per sample too
bcftools stats -s - Sample_X.filtered.vcf.gz > Sample_X.filtered.vcfstats.txt

# allele frequency spectrum, useful downstream for demography sanity checks
vcftools --gzvcf Sample_X.filtered.vcf.gz --freq2 --out Sample_X.filtered
