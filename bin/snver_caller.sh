#!/bin/bash

# SNVer variant caller script

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

# Create symbolic link to reference genome and index it
ln -sf "$REF_GENOME" reference.fasta
samtools faidx reference.fasta

# Run SNVer variant calling
snver -i "$BAM_FILE" -r reference.fasta -o "$BAM_BASENAME" -l "$COVERAGE_FILE" -bq "$MIN_BASE_QUALITY" -n "$PLOIDY"

# Process SNPs VCF
bcftools reheader -f reference.fasta.fai "${BAM_BASENAME}.filter.vcf" | \
    bcftools filter -e"QUAL<$MIN_SNP_QUAL" - | \
    bcftools annotate --force -x INFO,FORMAT - | bcftools view --min-alleles 2 --max-alleles 2 - | \
    bcftools norm --fasta-ref "$REF_GENOME" --atom-overlaps '.' --atomize -Oz -o "${BAM_BASENAME}.snps_snver.vcf.gz"

# Process indels VCF
bcftools reheader -f reference.fasta.fai "${BAM_BASENAME}.indel.filter.vcf" | \
    bcftools filter -e"QUAL<$MIN_SNP_QUAL" - | \
    bcftools annotate --force -x INFO,FORMAT - | bcftools view --min-alleles 2 --max-alleles 2 - | \
    bcftools norm --fasta-ref "$REF_GENOME" --atom-overlaps '.' --atomize -Oz -o "${BAM_BASENAME}.indels_snver.vcf.gz"

