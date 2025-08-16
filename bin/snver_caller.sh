#!/bin/bash

# SNVer variant caller script

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
    echo "[${timestamp}] [${level}] [SNVER_CALLER] [${BAM_BASENAME}] ${message}"
}

# Log process start
log_message "INFO" "Process started - SNVer variant calling"
log_message "INFO" "Parameters: BAM=${BAM_FILE}, REF=${REF_GENOME}, PLOIDY=${PLOIDY}, MIN_BQ=${MIN_BASE_QUALITY}, MIN_QUAL=${MIN_SNP_QUAL}"

# Record start time
START_TIME=$(date +%s)

# Create symbolic link to reference genome and index it
log_message "INFO" "Creating symbolic link to reference genome and indexing"
ln -sf "$REF_GENOME" reference.fasta
samtools faidx reference.fasta

if [ $? -eq 0 ]; then
    log_message "INFO" "Reference genome indexing completed successfully"
else
    log_message "ERROR" "Reference genome indexing failed"
    exit 1
fi

# Run SNVer variant calling
log_message "INFO" "Running SNVer variant calling"
snver -i "$BAM_FILE" -r reference.fasta -o "$BAM_BASENAME" -l "$COVERAGE_FILE" -bq "$MIN_BASE_QUALITY" -n "$PLOIDY"

if [ $? -eq 0 ]; then
    log_message "INFO" "SNVer variant calling completed successfully"
else
    log_message "ERROR" "SNVer variant calling failed"
    exit 1
fi

# Process SNPs VCF
log_message "INFO" "Processing SNPs VCF"
bcftools reheader -f reference.fasta.fai "${BAM_BASENAME}.filter.vcf" | \
    bcftools filter -e"QUAL<$MIN_SNP_QUAL" - | \
    bcftools annotate --force -x INFO,FORMAT - | bcftools view --min-alleles 2 --max-alleles 2 - | \
    bcftools norm --fasta-ref "$REF_GENOME" --atom-overlaps '.' --atomize -Oz -o "${BAM_BASENAME}.snps_snver.vcf.gz"

if [ $? -eq 0 ]; then
    log_message "INFO" "SNPs VCF processing completed successfully"
else
    log_message "ERROR" "SNPs VCF processing failed"
    exit 1
fi

# Process indels VCF
log_message "INFO" "Processing indels VCF"
bcftools reheader -f reference.fasta.fai "${BAM_BASENAME}.indel.filter.vcf" | \
    bcftools filter -e"QUAL<$MIN_SNP_QUAL" - | \
    bcftools annotate --force -x INFO,FORMAT - | bcftools view --min-alleles 2 --max-alleles 2 - | \
    bcftools norm --fasta-ref "$REF_GENOME" --atom-overlaps '.' --atomize -Oz -o "${BAM_BASENAME}.indels_snver.vcf.gz"

if [ $? -eq 0 ]; then
    log_message "INFO" "Indels VCF processing completed successfully"
else
    log_message "ERROR" "Indels VCF processing failed"
    exit 1
fi

# Calculate duration and log completion
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

log_message "INFO" "Process completed - SNPs and indels VCFs generated"
log_message "INFO" "Performance: ${DURATION} seconds"
log_message "INFO" "Output files: ${BAM_BASENAME}.snps_snver.vcf.gz, ${BAM_BASENAME}.indels_snver.vcf.gz"

# Clean up temporary files
log_message "INFO" "Cleaning up temporary files"
rm -f reference.fasta reference.fasta.fai "${BAM_BASENAME}.raw.vcf" "${BAM_BASENAME}.filter.vcf" "${BAM_BASENAME}.indel.filter.vcf" "${BAM_BASENAME}.indel.raw.vcf" "${BAM_BASENAME}.failed.log" 2>/dev/null || true
log_message "INFO" "Temporary files cleaned up"

