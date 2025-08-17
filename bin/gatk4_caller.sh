#!/bin/bash

# GATK4 variant caller script

set -e

# Input parameters
BAM_FILE="$1"
COVERAGE_FILE="$2"
REF_GENOME="$3"
PLOIDY="$4"
MIN_BASE_QUALITY="$5"
MIN_SNP_QUAL="$6"

# Validate required parameters
if [ $# -ne 6 ]; then
    echo "Usage: $0 <BAM_FILE> <COVERAGE_FILE> <REF_GENOME> <PLOIDY> <MIN_BASE_QUALITY> <MIN_SNP_QUAL>"
    echo "Error: Expected 6 parameters, got $#"
    exit 1
fi

# Validate that required files exist
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

# Extract base name from BAM file
BAM_BASENAME=$(basename "$BAM_FILE" .bam)

# Logging function
log_message() {
    local level="$1"
    local message="$2"
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
gatk HaplotypeCaller -I "$BAM_FILE" -R "$REF_GENOME" -mbq "$MIN_BASE_QUALITY" -O "${BAM_BASENAME}.gatk1.vcf.gz" -L "$COVERAGE_FILE" -ploidy "$PLOIDY" --do-not-run-physical-phasing true --smith-waterman FASTEST_AVAILABLE --create-output-variant-index false

if [ $? -eq 0 ]; then
    log_message "INFO" "GATK4 HaplotypeCaller completed successfully"
else
    log_message "ERROR" "GATK4 HaplotypeCaller failed"
    exit 1
fi

# Filter and process VCF
log_message "INFO" "Filtering and processing VCF with BCFtools"

# Step 1: Filter by quality
bcftools filter -Ou "${BAM_BASENAME}.gatk1.vcf.gz" -e"QUAL<$MIN_SNP_QUAL" > "${BAM_BASENAME}.gatk2.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools filter failed"
    exit 1
fi

# Step 2: Annotate and remove INFO/FORMAT fields
bcftools annotate -Ou --force -x INFO,FORMAT "${BAM_BASENAME}.gatk2.bcf" > "${BAM_BASENAME}.gatk3.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools annotate failed"
    exit 1
fi

# Step 3: Filter for biallelic variants
bcftools view -Ou --min-alleles 2 --max-alleles 2 "${BAM_BASENAME}.gatk3.bcf" > "${BAM_BASENAME}.gatk4.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools view (allele filter) failed"
    exit 1
fi

# Step 4: Normalize variants
bcftools norm -Ou --fasta-ref "$REF_GENOME" --atom-overlaps '.' --atomize "${BAM_BASENAME}.gatk4.bcf" > "${BAM_BASENAME}.gatk.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools norm failed"
    exit 1
fi

log_message "INFO" "VCF filtering and processing completed successfully"

# Extract SNPs and indels
log_message "INFO" "Extracting SNPs and indels from GATK4 output"
bcftools view -Ov -v snps "${BAM_BASENAME}.gatk.bcf" | bgzip > "${BAM_BASENAME}.snps_gatk.vcf.gz"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools view for SNPs failed"
    exit 1
fi

bcftools view -Ov -v indels "${BAM_BASENAME}.gatk.bcf" | bgzip > "${BAM_BASENAME}.indels_gatk.vcf.gz"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools view for indels failed"
    exit 1
fi

# Calculate duration and log completion
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

log_message "INFO" "Process completed - SNPs and indels extracted"
log_message "INFO" "Performance: ${DURATION} seconds"
log_message "INFO" "Output files: ${BAM_BASENAME}.snps_gatk.vcf.gz, ${BAM_BASENAME}.indels_gatk.vcf.gz"

