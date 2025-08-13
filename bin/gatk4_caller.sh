#!/bin/bash

# GATK4 variant caller script

set -e

# Input parameters
BAM_FILE=$1
COVERAGE_FILE=$2
REF_GENOME=$3
PLOIDY=$4
MIN_BASE_QUALITY=$5
MIN_SNP_QUAL=$6

# Extract base name from BAM file
BAM_BASENAME=$(basename "$BAM_FILE" .bam)

# Run GATK4 HaplotypeCaller with appropriate parameters based on ploidy
gatk HaplotypeCaller -I "$BAM_FILE" -R "$REF_GENOME" -mbq "$MIN_BASE_QUALITY" -O "${BAM_BASENAME}.gatk1.vcf.gz" -L "$COVERAGE_FILE" -ploidy "$PLOIDY" --do-not-run-physical-phasing true --smith-waterman FASTEST_AVAILABLE

bcftools filter "${BAM_BASENAME}.gatk1.vcf.gz" -e"QUAL<$MIN_SNP_QUAL" | \
    bcftools annotate --force -x INFO,FORMAT - | \
    bcftools view --min-alleles 2 --max-alleles 2 - | \
    bcftools norm --fasta-ref "$REF_GENOME" --atom-overlaps '.' --atomize -Oz -o "${BAM_BASENAME}.gatk.vcf.gz"

# Extract SNPs and indels
bcftools view -v snps -Oz -o "${BAM_BASENAME}.snps_gatk.vcf.gz" "${BAM_BASENAME}.gatk.vcf.gz"
bcftools view -v indels -Oz -o "${BAM_BASENAME}.indels_gatk.vcf.gz" "${BAM_BASENAME}.gatk.vcf.gz"

