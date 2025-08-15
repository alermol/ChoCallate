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

# Logging function
log_message() {
    local level=$1
    local message=$2
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
    var2vcf_valid.pl -q "$MIN_BASE_QUALITY" -N "$BAM_BASENAME" -E | \
    bcftools reheader -f "$REF_GENOME.fai" - | bcftools filter -e"QUAL<$MIN_SNP_QUAL" - | \
    bcftools annotate --force -x INFO,FORMAT - | bcftools view --min-alleles 2 --max-alleles 2 - | \
    bcftools norm --fasta-ref "$REF_GENOME" --atom-overlaps '.' --atomize -Oz -o "${BAM_BASENAME}.vardict.vcf.gz"

if [ $? -eq 0 ]; then
    log_message "INFO" "VarDict variant calling completed successfully"
else
    log_message "ERROR" "VarDict variant calling failed"
    exit 1
fi

# Extract SNPs and indels
log_message "INFO" "Extracting SNPs and indels from VarDict output"
bcftools view -v snps -Oz -o "${BAM_BASENAME}.snps_vardict.vcf.gz" "${BAM_BASENAME}.vardict.vcf.gz"
bcftools view -v indels -Oz -o "${BAM_BASENAME}.indels_vardict.vcf.gz" "${BAM_BASENAME}.vardict.vcf.gz"

# Calculate duration and log completion
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

log_message "INFO" "Process completed - SNPs and indels extracted"
log_message "INFO" "Performance: ${DURATION} seconds"
log_message "INFO" "Output files: ${BAM_BASENAME}.snps_vardict.vcf.gz, ${BAM_BASENAME}.indels_vardict.vcf.gz"

# Clean up temporary files
log_message "INFO" "Cleaning up temporary files"
rm -f "${BAM_BASENAME}.vardict.vcf.gz" 2>/dev/null || true
log_message "INFO" "Temporary files cleaned up"

