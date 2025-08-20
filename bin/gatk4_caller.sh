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
    echo "[${timestamp}] [${level}] [GATK4_CALLER] [${BAM_BASENAME}] ${message}"
}

log_message "INFO" "Process started - GATK4 HaplotypeCaller variant calling (ploidy=${PLOIDY}, min_qual=${MIN_SNP_QUAL})"

START_TIME=$(date +%s)

log_message "INFO" "Running GATK4 HaplotypeCaller"
gatk HaplotypeCaller -I "$BAM_FILE" -R "$REF_GENOME" -mbq "$MIN_BASE_QUALITY" -O "${BAM_BASENAME}.gatk1.vcf.gz" -L "$COVERAGE_FILE" -ploidy "$PLOIDY" --do-not-run-physical-phasing true --smith-waterman FASTEST_AVAILABLE --create-output-variant-index false

if [ $? -eq 0 ]; then
    log_message "INFO" "GATK4 HaplotypeCaller completed successfully"
else
    log_message "ERROR" "GATK4 HaplotypeCaller failed"
    exit 1
fi

log_message "INFO" "Indexing GATK4 VCF output with tabix"
tabix -p vcf -C "${BAM_BASENAME}.gatk1.vcf.gz"
if [ $? -ne 0 ]; then
    log_message "ERROR" "tabix indexing of GATK4 VCF failed"
    exit 1
fi

log_message "INFO" "Filtering and processing VCF with BCFtools"

bcftools filter -Ou "${BAM_BASENAME}.gatk1.vcf.gz" -e"QUAL<$MIN_SNP_QUAL" - | \
    bcftools annotate -Ou --force -x INFO,FORMAT - | \
    bcftools view -Ou --min-alleles 2 --max-alleles 2 - | \
    bcftools norm -Ob --fasta-ref "$REF_GENOME" --atom-overlaps '.' --atomize > "${BAM_BASENAME}.gatk.bcf"

if [ $? -ne 0 ]; then
    log_message "ERROR" "BCFtools filtering/processing pipeline failed"
    exit 1
fi

log_message "INFO" "VCF filtering and processing completed successfully"

log_message "INFO" "Extracting SNPs and indels from GATK4 output (compressed BCF)"
bcftools view -Ob -v snps "${BAM_BASENAME}.gatk.bcf" > "${BAM_BASENAME}.snps_gatk.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools view for SNPs failed"
    exit 1
fi

bcftools view -Ob -v indels "${BAM_BASENAME}.gatk.bcf" > "${BAM_BASENAME}.indels_gatk.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools view for indels failed"
    exit 1
fi

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

log_message "INFO" "Process completed - SNPs and indels extracted (${DURATION}s)"
log_message "INFO" "Output files: ${BAM_BASENAME}.snps_gatk.bcf, ${BAM_BASENAME}.indels_gatk.bcf"

