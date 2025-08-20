#!/bin/bash

set -e

BAM_FILE="$1"
COVERAGE_FILE="$2"
REF_GENOME="$3"
PLOIDY="$4"
MIN_BASE_QUALITY="$5"
MIN_SNP_QUAL="$6"

if [ $# -ne 6 ]; then
    echo "Usage: $0 <BAM_FILE> <COVERAGE_FILE> <REF_GENOME> <PLOIDY> <MIN_BASE_QUALITY> <MIN_SNP_QUAL>"
    echo "Error: Expected 6 parameters, got $#"
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
    echo "[${timestamp}] [${level}] [SNVER_CALLER] [${BAM_BASENAME}] ${message}"
}

log_message "INFO" "Process started - SNVer variant calling (ploidy=${PLOIDY}, min_qual=${MIN_SNP_QUAL})"

START_TIME=$(date +%s)

if [ $? -eq 0 ]; then
    log_message "INFO" "Reference genome indexing completed successfully"
else
    log_message "ERROR" "Reference genome indexing failed"
    exit 1
fi

log_message "INFO" "Running SNVer variant calling"
snver -i "$BAM_FILE" -r reference.fasta -o "$BAM_BASENAME" -l "$COVERAGE_FILE" -bq "$MIN_BASE_QUALITY" -n "$PLOIDY"

if [ $? -eq 0 ]; then
    log_message "INFO" "SNVer variant calling completed successfully"
else
    log_message "ERROR" "SNVer variant calling failed"
    exit 1
fi

log_message "INFO" "Processing SNPs VCF"

bcftools reheader -f reference.fasta.fai "${BAM_BASENAME}.filter.vcf" | \
    bcftools view -Ou - | \
    bcftools filter -Ou -e"QUAL<$MIN_SNP_QUAL" - | \
    bcftools annotate -Ou --force -x INFO,FORMAT - | \
    bcftools view -Ou --min-alleles 2 --max-alleles 2 - | \
    bcftools norm -Ob --fasta-ref reference.fasta --atom-overlaps '.' --atomize - > "${BAM_BASENAME}.snps.bcf"

if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools SNPs processing pipeline failed"
    exit 1
fi

log_message "INFO" "SNPs VCF processing completed successfully"

log_message "INFO" "Processing indels VCF"

bcftools reheader -f reference.fasta.fai "${BAM_BASENAME}.indel.filter.vcf" | \
    bcftools view -Ou - | \
    bcftools filter -Ou -e"QUAL<$MIN_SNP_QUAL" - | \
    bcftools annotate -Ou --force -x INFO,FORMAT - | \
    bcftools view -Ou --min-alleles 2 --max-alleles 2 - | \
    bcftools norm -Ob --fasta-ref reference.fasta --atom-overlaps '.' --atomize - > "${BAM_BASENAME}.indels.bcf"

if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools indels processing pipeline failed"
    exit 1
fi

log_message "INFO" "Indels VCF processing completed successfully"

log_message "INFO" "Generating final compressed BCF output files"
bcftools view -Ob "${BAM_BASENAME}.snps.bcf" > "${BAM_BASENAME}.snps_snver.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "Final SNPs output generation failed"
    exit 1
fi

bcftools view -Ob "${BAM_BASENAME}.indels.bcf" > "${BAM_BASENAME}.indels_snver.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "Final indels output generation failed"
    exit 1
fi

log_message "INFO" "Final compressed BCF output files generated successfully"

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

log_message "INFO" "Process completed - SNPs and indels BCFs generated (${DURATION}s)"
log_message "INFO" "Output files: ${BAM_BASENAME}.snps_snver.bcf, ${BAM_BASENAME}.indels_snver.bcf"