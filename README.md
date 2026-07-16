# Variant-calling for population genomics 

- Multi-sample mapping and variant calling against a de novo reference, with population-level
filtering.
- This step depends on what type and quality of variants you want. I worked mostly on population
genomics with low/ medium coverage per sample. For that reason, I prefer to filter SNPs of high
quality. However, if you work with mutations, filtering indels can be interesting.

## Pipeline

- `bwa-mem2` index of the reference, mapping of each sample's trimmed paired reads
- `samtools sort` / `markdup` / `index`, `flagstat` for mapping rate sanity checks per sample
- `bcftools mpileup` + `call` for joint calling across all samples
- filtering: biallelic sites only, SNPs only, present in at least 80% of samples, QUAL >= 30
- `bcftools stats -s` for statistic assessment - Sample_X.filtered.vcf.gz > Sample_X.filtered.vcfstats.txt

## Input

- `hifiasm_Sample_X.p_ctg.fasta`: reference assembly
- `Sample_*_R1.paired.fastq.gz`, `Sample_*_R2.paired.fastq.gz`: trimmed reads, from
  `02-qc-trimming`, one pair per sample in the `SAMPLES` array

## Output

- `Sample_*.dedup.bam`, `Sample_*.flagstat.txt`: per-sample mapping, for mapping-rate QC
- `Sample_X.raw.vcf.gz`: joint-called, unfiltered variants
- `Sample_X.filtered.vcf.gz` (+ `.tbi`): biallelic SNPs, QUAL >= 30, max 20% missingness
- `Sample_X.filtered.vcfstats.txt`: variants statistics

`Sample_X.filtered.vcf.gz` and the per-sample `.dedup.bam` files are the inputs reused in
`04-reference-consensus` and `05-psmc-demographic-history`.

## Usage

```
bwa-mem2 mem -t 32 -R "@RG\tID:Sample_X\tSM:Sample_X\tPL:ILLUMINA" \
    hifiasm_reference.p_ctg.fasta Sample_X_R1.paired.fastq.gz Sample_X_R2.paired.fastq.gz
samtools markdup -@ 32 Sample_X.sorted.bam Sample_X.dedup.bam
samtools index Sample_X.dedup.bam
samtools flagstat Sample_X.dedup.bam > Sample_X.flagstat.txt

bcftools mpileup -f hifiasm_reference.p_ctg.fasta -b Sample_X.flagstat.txt -a AD,DP --threads 32 -Ou \
  | bcftools call -mv --threads 32 -Oz -o Sample_X.raw.vcf.gz
tabix -p vcf Sample_X.raw.vcf.gz

bcftools view -m2 -M2 -v snps Sample_X.raw.vcf.gz -Oz \
  | bcftools view -e 'F_MISSING > 0.2' -Oz \                        #in at least 80% of samples
  | bcftools filter -e 'QUAL<30' -Oz -o Sample_X.filtered.vcf.gz    #QUAL >= 30
tabix -p vcf Sample_X.filtered.vcf.gz

bcftools stats -s - Sample_X.filtered.vcf.gz > Sample_X.filtered.vcfstats.txt
```

## Notes

I prefer bcftools over GATK for non-model organisms, since there's no known-sites VCF available
for BQSR and it's considerably lighter to run on a shared HPC. If you go with GATK instead you'll
need a sequence dictionary (`picard CreateSequenceDictionary`) and `HaplotypeCaller` run per
sample before joint genotyping with `GenomicsDBImport` / `GenotypeGVCFs`.
For more details, check  `broadinstitute/gatk`
