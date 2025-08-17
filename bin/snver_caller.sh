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

log_message "INFO" "Creating symbolic link to reference genome and indexing"
ln -sf "$REF_GENOME" reference.fasta
samtools faidx reference.fasta

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

bcftools reheader -f reference.fasta.fai "${BAM_BASENAME}.filter.vcf" | bgzip > "${BAM_BASENAME}.snver1.vcf.gz"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools reheader for SNPs failed"
    exit 1
fi

bcftools filter -Ou -e"QUAL<$MIN_SNP_QUAL" "${BAM_BASENAME}.snver1.vcf.gz" > "${BAM_BASENAME}.snver2.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools filter for SNPs failed"
    exit 1
fi

bcftools annotate -Ou --force -x INFO,FORMAT "${BAM_BASENAME}.snver2.bcf" > "${BAM_BASENAME}.snver3.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools annotate for SNPs failed"
    exit 1
fi

bcftools view -Ou --min-alleles 2 --max-alleles 2 "${BAM_BASENAME}.snver3.bcf" > "${BAM_BASENAME}.snver4.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools view (allele filter) for SNPs failed"
    exit 1
fi

bcftools norm -Ou --fasta-ref reference.fasta --atom-overlaps '.' --atomize "${BAM_BASENAME}.snver4.bcf" > "${BAM_BASENAME}.snps_snver.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools norm for SNPs failed"
    exit 1
fi

log_message "INFO" "SNPs VCF processing completed successfully"

log_message "INFO" "Processing indels VCF"

bcftools reheader -f reference.fasta.fai "${BAM_BASENAME}.indel.filter.vcf" | bgzip > "${BAM_BASENAME}.snver1.vcf.gz"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools reheader for indels failed"
    exit 1
fi

bcftools filter -Ou -e"QUAL<$MIN_SNP_QUAL" "${BAM_BASENAME}.snver1.vcf.gz" > "${BAM_BASENAME}.snver2.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools filter for indels failed"
    exit 1
fi

bcftools annotate -Ou --force -x INFO,FORMAT "${BAM_BASENAME}.snver2.bcf" > "${BAM_BASENAME}.snver3.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools annotate for indels failed"
    exit 1
fi

bcftools view -Ou --min-alleles 2 --max-alleles 2 "${BAM_BASENAME}.snver3.bcf" > "${BAM_BASENAME}.snver4.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools view (allele filter) for indels failed"
    exit 1
fi

bcftools norm -Ou --fasta-ref reference.fasta --atom-overlaps '.' --atomize "${BAM_BASENAME}.snver4.bcf" > "${BAM_BASENAME}.indels_snver.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools norm for indels failed"
    exit 1
fi

log_message "INFO" "Indels VCF processing completed successfully"

log_message "INFO" "Generating final compressed output files"
bcftools view -Ov "${BAM_BASENAME}.snps_snver.bcf" | bgzip > "${BAM_BASENAME}.snps_snver.vcf.gz"
if [ $? -ne 0 ]; then
    log_message "ERROR" "Final SNPs output generation failed"
    exit 1
fi

bcftools view -Ov "${BAM_BASENAME}.indels_snver.bcf" | bgzip > "${BAM_BASENAME}.indels_snver.vcf.gz"
if [ $? -ne 0 ]; then
    log_message "ERROR" "Final indels output generation failed"
    exit 1
fi

log_message "INFO" "Final compressed output files generated successfully"

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

log_message "INFO" "Process completed - SNPs and indels VCFs generated (${DURATION}s)"
log_message "INFO" "Output files: ${BAM_BASENAME}.snps_snver.vcf.gz, ${BAM_BASENAME}.indels_snver.vcf.gz"