#!/bin/bash

set -e

BAM_FILE="$1"
COVERAGE_FILE="$2"
REF_GENOME="$3"
PLOIDY="$4"
MIN_BASE_QUALITY="$5"
MIN_SNP_QUAL="$6"
THREADS="$7"
MPILEUP_EXTRA_ARGS="$8"
CALL_EXTRA_ARGS="$9"

if [ $# -lt 7 ] || [ $# -gt 9 ]; then
    echo "Usage: $0 <BAM_FILE> <COVERAGE_FILE> <REF_GENOME> <PLOIDY> <MIN_BASE_QUALITY> <MIN_SNP_QUAL> <THREADS> <MPILEUP_EXTRA_ARGS> <CALL_EXTRA_ARGS>"
    echo "Error: Expected 9 parameters, got $#"
    exit 1
fi

if [ ! -f "$BAM_FILE" ]; then
    echo "Error: BAM file '$BAM_FILE' not found"
    exit 1
fi

if [ ! -f "$COVERAGE_FILE" ]; then
    echo "Error: Coverage file '$COVERAGE_FILE' not found"
    exit 1
fi

if [ ! -f "$REF_GENOME" ]; then
    echo "Error: Reference genome '$REF_GENOME' not found"
    exit 1
fi

BAM_BASENAME=$(basename "$BAM_FILE" .bam)

log_message() {
    local level="$1"
    local message="$2"
    local timestamp=$(date -Iseconds)
    echo "[${timestamp}] [${level}] [BCFTOOLS_CALLER] [${BAM_BASENAME}] ${message}"
}

log_message "INFO" "Process started - BCFtools variant calling (ploidy=${PLOIDY}, min_qual=${MIN_SNP_QUAL}, threads=${THREADS})"

START_TIME=$(date +%s)

log_message "INFO" "Running BCFtools mpileup and variant calling"
bcftools mpileup -Ou --count-orphans --fasta-ref "$REF_GENOME" --threads "$THREADS" --max-depth 250 \
    --min-BQ "$MIN_BASE_QUALITY" --regions-file "$COVERAGE_FILE" ${MPILEUP_EXTRA_ARGS} "$BAM_FILE" | \
    bcftools call -Ou --multiallelic-caller --threads "$THREADS" ${CALL_EXTRA_ARGS} | \
    bcftools filter -Ou -e"QUAL<$MIN_SNP_QUAL" - | \
    bcftools annotate -Ou --force -x INFO,FORMAT - | \
    bcftools view -Ou --min-alleles 2 --max-alleles 2 - | \
    bcftools norm -Ob --fasta-ref "$REF_GENOME" --atom-overlaps '.' --atomize > "${BAM_BASENAME}.bcftools.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools variant calling pipeline failed"
    exit 1
fi

log_message "INFO" "BCFtools variant calling completed successfully"

log_message "INFO" "Extracting SNPs and indels from BCFtools output"
bcftools view -Ob -v snps "${BAM_BASENAME}.bcftools.bcf" > "${BAM_BASENAME}.snps_bcftools.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools view for SNPs failed"
    exit 1
fi

bcftools view -Ob -v indels "${BAM_BASENAME}.bcftools.bcf" > "${BAM_BASENAME}.indels_bcftools.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools view for indels failed"
    exit 1
fi

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

log_message "INFO" "Process completed - SNPs and indels extracted (${DURATION}s)"
log_message "INFO" "Output files: ${BAM_BASENAME}.snps_bcftools.bcf, ${BAM_BASENAME}.indels_bcftools.bcf"


