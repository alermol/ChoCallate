#!/bin/bash

set -e

BAM_FILE="$1"
COVERAGE_FILE="$2"
REF_GENOME="$3"
PLOIDY="$4"
MIN_BASE_QUALITY="$5"
MIN_SNP_QUAL="$6"
THREADS="$7"

if [ $# -ne 7 ]; then
    echo "Usage: $0 <BAM_FILE> <COVERAGE_FILE> <REF_GENOME> <PLOIDY> <MIN_BASE_QUALITY> <MIN_SNP_QUAL> <THREADS>"
    echo "Error: Expected 7 parameters, got $#"
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
    echo "[${timestamp}] [${level}] [VARDICT_CALLER] [${BAM_BASENAME}] ${message}"
}

log_message "INFO" "Process started - VarDict variant calling (ploidy=${PLOIDY}, min_qual=${MIN_SNP_QUAL}, threads=${THREADS})"

START_TIME=$(date +%s)

log_message "INFO" "Running VarDict variant calling"
vardict-java -G "$REF_GENOME" -N "$BAM_BASENAME" -b "$BAM_FILE" -fisher -th "$THREADS" \
    -VS SILENT --nosv -k 0 -q "$MIN_BASE_QUALITY" -c 1 -S 2 -E 3 -g 4 "$COVERAGE_FILE" | \
    var2vcf_valid.pl -S -q "$MIN_BASE_QUALITY" -N "$BAM_BASENAME" -E | \
    bcftools reheader -f "$REF_GENOME.fai" | \
    bcftools view -Ou - | \
    bcftools filter -Ou -e"QUAL<$MIN_SNP_QUAL" - | \
    bcftools annotate -Ou --force -x INFO,FORMAT - | \
    bcftools view -Ou --min-alleles 2 --max-alleles 2 - | \
    bcftools norm -Ob --fasta-ref "$REF_GENOME" --atom-overlaps '.' --atomize - > "${BAM_BASENAME}.vardict.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "VarDict variant calling pipeline failed"
    exit 1
fi

log_message "INFO" "VarDict variant calling completed successfully"

log_message "INFO" "Extracting SNPs and indels from VarDict output"
bcftools view -Ob -v snps "${BAM_BASENAME}.vardict.bcf" > "${BAM_BASENAME}.snps_vardict.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools view for SNPs failed"
    exit 1
fi

bcftools view -Ob -v indels "${BAM_BASENAME}.vardict.bcf" > "${BAM_BASENAME}.indels_vardict.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools view for indels failed"
    exit 1
fi

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

log_message "INFO" "Process completed - SNPs and indels extracted (${DURATION}s)"
log_message "INFO" "Output files: ${BAM_BASENAME}.snps_vardict.bcf, ${BAM_BASENAME}.indels_vardict.bcf"