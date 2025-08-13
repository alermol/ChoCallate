#!/bin/bash

# BCFtools variant caller script

set -e

# Input parameters
BAM_FILE=$1
COVERAGE_FILE=$2
REF_GENOME=$3
PLOIDY=$4
MIN_BASE_QUALITY=$5
MIN_SNP_QUAL=$6
THREADS=$7

# Extract base name from BAM file
BAM_BASENAME=$(basename "$BAM_FILE" .bam)

# Run BCFtools mpileup and call variants
bcftools mpileup -Ou --count-orphans --fasta-ref "$REF_GENOME" --threads "$THREADS" --max-depth 250 \
    --min-BQ "$MIN_BASE_QUALITY" --regions-file "$COVERAGE_FILE" "$BAM_FILE" | \
    bcftools call -Ov --multiallelic-caller --threads "$THREADS" | \
    bcftools filter -e"QUAL<$MIN_SNP_QUAL" - | \
    bcftools annotate --force -x INFO,FORMAT - | bcftools view --min-alleles 2 --max-alleles 2 - | \
    bcftools norm --fasta-ref "$REF_GENOME" --atom-overlaps '.' --atomize -Oz -o "${BAM_BASENAME}.bcftools"

# Extract SNPs and indels
bcftools view -v snps -Oz -o "${BAM_BASENAME}.snps_bcftools.vcf.gz" "${BAM_BASENAME}.bcftools"
bcftools view -v indels -Oz -o "${BAM_BASENAME}.indels_bcftools.vcf.gz" "${BAM_BASENAME}.bcftools"


