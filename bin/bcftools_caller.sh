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

# Logging function
log_message() {
    local level=$1
    local message=$2
    local timestamp=$(date -Iseconds)
    echo "[${timestamp}] [${level}] [BCFTOOLS_CALLER] [${BAM_BASENAME}] ${message}"
}

# Log process start
log_message "INFO" "Process started - BCFtools variant calling"
log_message "INFO" "Parameters: BAM=${BAM_FILE}, REF=${REF_GENOME}, PLOIDY=${PLOIDY}, MIN_BQ=${MIN_BASE_QUALITY}, MIN_QUAL=${MIN_SNP_QUAL}, THREADS=${THREADS}"

# Record start time
START_TIME=$(date +%s)

# Run BCFtools mpileup and call variants
log_message "INFO" "Running BCFtools mpileup and variant calling"
bcftools mpileup -Ou --count-orphans --fasta-ref "$REF_GENOME" --threads "$THREADS" --max-depth 250 \
    --min-BQ "$MIN_BASE_QUALITY" --regions-file "$COVERAGE_FILE" "$BAM_FILE" | \
    bcftools call -Ov --multiallelic-caller --threads "$THREADS" | \
    bcftools filter -e"QUAL<$MIN_SNP_QUAL" - | \
    bcftools annotate --force -x INFO,FORMAT - | bcftools view --min-alleles 2 --max-alleles 2 - | \
    bcftools norm --fasta-ref "$REF_GENOME" --atom-overlaps '.' --atomize -Oz -o "${BAM_BASENAME}.bcftools"

if [ $? -eq 0 ]; then
    log_message "INFO" "BCFtools variant calling completed successfully"
else
    log_message "ERROR" "BCFtools variant calling failed"
    exit 1
fi

# Extract SNPs and indels
log_message "INFO" "Extracting SNPs and indels from BCFtools output"
bcftools view -v snps -Oz -o "${BAM_BASENAME}.snps_bcftools.vcf.gz" "${BAM_BASENAME}.bcftools"
bcftools view -v indels -Oz -o "${BAM_BASENAME}.indels_bcftools.vcf.gz" "${BAM_BASENAME}.bcftools"

# Calculate duration and log completion
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

log_message "INFO" "Process completed - SNPs and indels extracted"
log_message "INFO" "Performance: ${DURATION} seconds"
log_message "INFO" "Output files: ${BAM_BASENAME}.snps_bcftools.vcf.gz, ${BAM_BASENAME}.indels_bcftools.vcf.gz"


