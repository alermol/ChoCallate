#!/bin/bash

# BCFtools variant caller script

set -e

# Input parameters
BAM_FILE="$1"
COVERAGE_FILE="$2"
REF_GENOME="$3"
PLOIDY="$4"
MIN_BASE_QUALITY="$5"
MIN_SNP_QUAL="$6"
THREADS="$7"

# Validate required parameters
if [ $# -ne 7 ]; then
    echo "Usage: $0 <BAM_FILE> <COVERAGE_FILE> <REF_GENOME> <PLOIDY> <MIN_BASE_QUALITY> <MIN_SNP_QUAL> <THREADS>"
    echo "Error: Expected 7 parameters, got $#"
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
    bcftools call -Ou --multiallelic-caller --threads "$THREADS" > "${BAM_BASENAME}.bcftools1.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools mpileup and call failed"
    exit 1
fi

# Step 3: Filter by quality
bcftools filter -Ou -e"QUAL<$MIN_SNP_QUAL" "${BAM_BASENAME}.bcftools1.bcf" > "${BAM_BASENAME}.bcftools2.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools filter failed"
    exit 1
fi

# Step 4: Annotate and remove INFO/FORMAT fields
bcftools annotate -Ou --force -x INFO,FORMAT "${BAM_BASENAME}.bcftools2.bcf" > "${BAM_BASENAME}.bcftools3.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools annotate failed"
    exit 1
fi

# Step 5: Filter for biallelic variants
bcftools view -Ou --min-alleles 2 --max-alleles 2 "${BAM_BASENAME}.bcftools3.bcf" > "${BAM_BASENAME}.bcftools4.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools view (allele filter) failed"
    exit 1
fi

# Step 6: Normalize variants
bcftools norm -Ou --fasta-ref "$REF_GENOME" --atom-overlaps '.' --atomize "${BAM_BASENAME}.bcftools4.bcf" > "${BAM_BASENAME}.bcftools.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools norm failed"
    exit 1
fi

log_message "INFO" "BCFtools variant calling completed successfully"

# Extract SNPs and indels
log_message "INFO" "Extracting SNPs and indels from BCFtools output"
bcftools view -Ov -v snps "${BAM_BASENAME}.bcftools.bcf" | bgzip > "${BAM_BASENAME}.snps_bcftools.vcf.gz"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools view for SNPs failed"
    exit 1
fi

bcftools view -Ov -v indels "${BAM_BASENAME}.bcftools.bcf" | bgzip > "${BAM_BASENAME}.indels_bcftools.vcf.gz"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools view for indels failed"
    exit 1
fi

# Calculate duration and log completion
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

log_message "INFO" "Process completed - SNPs and indels extracted"
log_message "INFO" "Performance: ${DURATION} seconds"
log_message "INFO" "Output files: ${BAM_BASENAME}.snps_bcftools.vcf.gz, ${BAM_BASENAME}.indels_bcftools.vcf.gz"


