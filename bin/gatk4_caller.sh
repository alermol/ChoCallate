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

# Logging function
log_message() {
    local level=$1
    local message=$2
    local timestamp=$(date -Iseconds)
    echo "[${timestamp}] [${level}] [GATK4_CALLER] [${BAM_BASENAME}] ${message}"
}

# Log process start
log_message "INFO" "Process started - GATK4 HaplotypeCaller variant calling"
log_message "INFO" "Parameters: BAM=${BAM_FILE}, REF=${REF_GENOME}, PLOIDY=${PLOIDY}, MIN_BQ=${MIN_BASE_QUALITY}, MIN_QUAL=${MIN_SNP_QUAL}"

# Record start time
START_TIME=$(date +%s)

# Run GATK4 HaplotypeCaller with appropriate parameters based on ploidy
log_message "INFO" "Running GATK4 HaplotypeCaller"
gatk HaplotypeCaller -I "$BAM_FILE" -R "$REF_GENOME" -mbq "$MIN_BASE_QUALITY" -O "${BAM_BASENAME}.gatk1.vcf.gz" -L "$COVERAGE_FILE" -ploidy "$PLOIDY" --do-not-run-physical-phasing true --smith-waterman FASTEST_AVAILABLE

if [ $? -eq 0 ]; then
    log_message "INFO" "GATK4 HaplotypeCaller completed successfully"
else
    log_message "ERROR" "GATK4 HaplotypeCaller failed"
    exit 1
fi

# Filter and process VCF
log_message "INFO" "Filtering and processing VCF with BCFtools"
bcftools filter "${BAM_BASENAME}.gatk1.vcf.gz" -e"QUAL<$MIN_SNP_QUAL" | \
    bcftools annotate --force -x INFO,FORMAT - | \
    bcftools view --min-alleles 2 --max-alleles 2 - | \
    bcftools norm --fasta-ref "$REF_GENOME" --atom-overlaps '.' --atomize -Oz -o "${BAM_BASENAME}.gatk.vcf.gz"

if [ $? -eq 0 ]; then
    log_message "INFO" "VCF filtering and processing completed successfully"
else
    log_message "ERROR" "VCF filtering and processing failed"
    exit 1
fi

# Extract SNPs and indels
log_message "INFO" "Extracting SNPs and indels from GATK4 output"
bcftools view -v snps -Oz -o "${BAM_BASENAME}.snps_gatk.vcf.gz" "${BAM_BASENAME}.gatk.vcf.gz"
bcftools view -v indels -Oz -o "${BAM_BASENAME}.indels_gatk.vcf.gz" "${BAM_BASENAME}.gatk.vcf.gz"

# Calculate duration and log completion
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

log_message "INFO" "Process completed - SNPs and indels extracted"
log_message "INFO" "Performance: ${DURATION} seconds"
log_message "INFO" "Output files: ${BAM_BASENAME}.snps_gatk.vcf.gz, ${BAM_BASENAME}.indels_gatk.vcf.gz"

