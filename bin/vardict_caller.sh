#!/bin/bash

# VarDict variant caller script

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

# Run VarDict variant calling
vardict-java -G "$REF_GENOME" -N "$BAM_BASENAME" -b "$BAM_FILE" -fisher -th "$THREADS" \
    -VS SILENT --nosv -k 0 -q "$MIN_BASE_QUALITY" -c 1 -S 2 -E 3 -g 4 "$COVERAGE_FILE" | \
    var2vcf_valid.pl -q "$MIN_BASE_QUALITY" -N "$BAM_BASENAME" -E | \
    bcftools reheader -f "$REF_GENOME.fai" - | bcftools filter -e"QUAL<$MIN_SNP_QUAL" - | \
    bcftools annotate --force -x INFO,FORMAT - | bcftools view --min-alleles 2 --max-alleles 2 - | \
    bcftools norm --fasta-ref "$REF_GENOME" --atom-overlaps '.' --atomize -Oz -o "${BAM_BASENAME}.vardict.vcf.gz"

# Extract SNPs and indels
bcftools view -v snps -Oz -o "${BAM_BASENAME}.snps_vardict.vcf.gz" "${BAM_BASENAME}.vardict.vcf.gz"
bcftools view -v indels -Oz -o "${BAM_BASENAME}.indels_vardict.vcf.gz" "${BAM_BASENAME}.vardict.vcf.gz"

