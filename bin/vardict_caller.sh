#!/bin/bash

# VarDict variant caller script

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
    echo "[${timestamp}] [${level}] [VARDICT_CALLER] [${BAM_BASENAME}] ${message}"
}

# Log process start
log_message "INFO" "Process started - VarDict variant calling"
log_message "INFO" "Parameters: BAM=${BAM_FILE}, REF=${REF_GENOME}, PLOIDY=${PLOIDY}, MIN_BQ=${MIN_BASE_QUALITY}, MIN_QUAL=${MIN_SNP_QUAL}, THREADS=${THREADS}"

# Record start time
START_TIME=$(date +%s)

# Run VarDict variant calling
log_message "INFO" "Running VarDict variant calling"
vardict-java -G "$REF_GENOME" -N "$BAM_BASENAME" -b "$BAM_FILE" -fisher -th "$THREADS" \
    -VS SILENT --nosv -k 0 -q "$MIN_BASE_QUALITY" -c 1 -S 2 -E 3 -g 4 "$COVERAGE_FILE" | \
    var2vcf_valid.pl -S -q "$MIN_BASE_QUALITY" -N "$BAM_BASENAME" -E > "${BAM_BASENAME}.vardict1.vcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "VarDict variant calling failed"
    exit 1
fi

# Split remaining pipes into individual steps
bcftools reheader -f "$REF_GENOME.fai" "${BAM_BASENAME}.vardict1.vcf" > "${BAM_BASENAME}.vardict2.vcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools reheader failed"
    exit 1
fi

bcftools filter -Ou -e"QUAL<$MIN_SNP_QUAL" "${BAM_BASENAME}.vardict2.vcf" > "${BAM_BASENAME}.vardict3.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools filter failed"
    exit 1
fi

bcftools annotate -Ou --force -x INFO,FORMAT "${BAM_BASENAME}.vardict3.bcf" > "${BAM_BASENAME}.vardict4.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools annotate failed"
    exit 1
fi

bcftools view -Ou --min-alleles 2 --max-alleles 2 "${BAM_BASENAME}.vardict4.bcf" > "${BAM_BASENAME}.vardict5.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools view (allele filter) failed"
    exit 1
fi

bcftools norm -Ou --fasta-ref "$REF_GENOME" --atom-overlaps '.' --atomize "${BAM_BASENAME}.vardict5.bcf" > "${BAM_BASENAME}.vardict.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools norm failed"
    exit 1
fi

log_message "INFO" "VarDict variant calling completed successfully"

# Extract SNPs and indels
log_message "INFO" "Extracting SNPs and indels from VarDict output"
bcftools view -Ov -v snps "${BAM_BASENAME}.vardict.bcf" | bgzip > "${BAM_BASENAME}.snps_vardict.vcf.gz"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools view for SNPs failed"
    exit 1
fi

bcftools view -Ov -v indels "${BAM_BASENAME}.vardict.bcf" | bgzip > "${BAM_BASENAME}.indels_vardict.vcf.gz"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools view for indels failed"
    exit 1
fi

# Calculate duration and log completion
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

log_message "INFO" "Process completed - SNPs and indels extracted"
log_message "INFO" "Performance: ${DURATION} seconds"
log_message "INFO" "Output files: ${BAM_BASENAME}.snps_vardict.vcf.gz, ${BAM_BASENAME}.indels_vardict.vcf.gz"