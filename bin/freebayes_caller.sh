#!/bin/bash

# FreeBayes variant caller script

set -e

# Input parameters
BAM_FILE=$1
COVERAGE_FILE=$2
REF_GENOME=$3
PLOIDY=$4
READS_SOURCE=$5
MIN_BASE_QUALITY=$6
MIN_SNP_QUAL=$7

# Extract base name from BAM file
BAM_BASENAME=$(basename "$BAM_FILE" .bam)

# Run FreeBayes with appropriate parameters based on reads source
if [[ "$READS_SOURCE" == "gbs" ]]; then
    freebayes --fasta-reference "$REF_GENOME" --targets "$COVERAGE_FILE" --dont-left-align-indels --ploidy "$PLOIDY" \
        --use-best-n-alleles 4 --min-alternate-qsum "$MIN_BASE_QUALITY" --hwe-priors-off --no-population-priors \
        --binomial-obs-priors-off --allele-balance-priors-off --min-base-quality "$MIN_BASE_QUALITY" \
        --haplotype-length -1 --throw-away-complex-obs --no-partial-observations --bam "$BAM_FILE" --limit-coverage 250 | \
        bcftools filter -e"QUAL<$MIN_SNP_QUAL" - | \
        bcftools view --min-alleles 2 --max-alleles 2 - | bcftools annotate --force -x INFO,FORMAT - | \
        bcftools norm --fasta-ref "$REF_GENOME" --atom-overlaps '.' --atomize -Oz -o "${BAM_BASENAME}.freebayes.vcf.gz"
else
    freebayes --fasta-reference "$REF_GENOME" --targets "$COVERAGE_FILE" --dont-left-align-indels --ploidy "$PLOIDY" \
        --use-best-n-alleles 4 --min-alternate-qsum "$MIN_BASE_QUALITY" --hwe-priors-off --no-population-priors \
        --allele-balance-priors-off --min-base-quality "$MIN_BASE_QUALITY" \
        --haplotype-length -1 --throw-away-complex-obs --no-partial-observations --bam "$BAM_FILE" --limit-coverage 250 | \
        bcftools filter -e"QUAL<$MIN_SNP_QUAL" - | \
        bcftools view --min-alleles 2 --max-alleles 2 - | bcftools annotate --force -x INFO,FORMAT - | \
        bcftools norm --fasta-ref "$REF_GENOME" --atom-overlaps '.' --atomize -Oz -o "${BAM_BASENAME}.freebayes.vcf.gz"
fi

# Extract SNPs and indels
bcftools view -v snps -Oz -o "${BAM_BASENAME}.snps_freebayes.vcf.gz" "${BAM_BASENAME}.freebayes.vcf.gz"
bcftools view -v indels -Oz -o "${BAM_BASENAME}.indels_freebayes.vcf.gz" "${BAM_BASENAME}.freebayes.vcf.gz"

