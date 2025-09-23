#!/bin/bash

set -e

BAM_FILE="$1"
COVERAGE_FILE="$2"
REF_GENOME="$3"
PLOIDY="$4"
READS_SOURCE="$5"
MIN_BASE_QUALITY="$6"
MIN_SNP_QUAL="$7"
FREEBAYES_EXTRA_ARGS="$8"

if [ $# -lt 7 ] || [ $# -gt 8 ]; then
    echo "Usage: $0 <BAM_FILE> <COVERAGE_FILE> <REF_GENOME> <PLOIDY> <READS_SOURCE> <MIN_BASE_QUALITY> <MIN_SNP_QUAL> <FREEBAYES_EXTRA_ARGS>"
    echo "Error: Expected 8 parameters, got $#"
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
    echo "[${timestamp}] [${level}] [FREEBAYES_CALLER] [${BAM_BASENAME}] ${message}"
}

log_message "INFO" "Process started - FreeBayes variant calling (ploidy=${PLOIDY}, reads_source=${READS_SOURCE}, min_qual=${MIN_SNP_QUAL})"

START_TIME=$(date +%s)

log_message "INFO" "Running FreeBayes variant calling with ${READS_SOURCE} parameters"
if [[ "$READS_SOURCE" == "gbs" ]]; then
    freebayes --fasta-reference "$REF_GENOME" --targets "$COVERAGE_FILE" --dont-left-align-indels --ploidy "$PLOIDY" \
        --use-best-n-alleles 4 --min-alternate-qsum "$MIN_BASE_QUALITY" --hwe-priors-off --no-population-priors \
        --binomial-obs-priors-off --allele-balance-priors-off --min-base-quality "$MIN_BASE_QUALITY" ${FREEBAYES_EXTRA_ARGS} \
        --haplotype-length -1 --throw-away-complex-obs --no-partial-observations --bam "$BAM_FILE" --limit-coverage 250 | \
        bcftools filter -e"QUAL<${MIN_SNP_QUAL}" -Ou | \
        bcftools view -Ou --min-alleles 2 --max-alleles 2 | \
        bcftools annotate -Ou --force -x INFO,FORMAT | \
        bcftools norm -Ob --fasta-ref "$REF_GENOME" --atom-overlaps '.' --atomize -o "${BAM_BASENAME}.freebayes.bcf"
    if [ $? -ne 0 ]; then
        log_message "ERROR" "FreeBayes variant calling and processing failed"
        exit 1
    fi
else
    freebayes --fasta-reference "$REF_GENOME" --targets "$COVERAGE_FILE" --dont-left-align-indels --ploidy "$PLOIDY" \
        --use-best-n-alleles 4 --min-alternate-qsum "$MIN_BASE_QUALITY" --hwe-priors-off --no-population-priors \
        --allele-balance-priors-off --min-base-quality "$MIN_BASE_QUALITY" ${FREEBAYES_EXTRA_ARGS} \
        --haplotype-length -1 --throw-away-complex-obs --no-partial-observations --bam "$BAM_FILE" --limit-coverage 250 | \
        bcftools filter -e"QUAL<${MIN_SNP_QUAL}" -Ou | \
        bcftools view -Ou --min-alleles 2 --max-alleles 2 | \
        bcftools annotate -Ou --force -x INFO,FORMAT | \
        bcftools norm -Ob --fasta-ref "$REF_GENOME" --atom-overlaps '.' --atomize -o "${BAM_BASENAME}.freebayes.bcf"
    if [ $? -ne 0 ]; then
        log_message "ERROR" "FreeBayes variant calling and processing failed"
        exit 1
    fi
fi

log_message "INFO" "FreeBayes variant calling completed successfully"

log_message "INFO" "Extracting SNPs and indels from FreeBayes output"
bcftools view -Ob -v snps "${BAM_BASENAME}.freebayes.bcf" -o "${BAM_BASENAME}.snps_freebayes.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools view for SNPs failed"
    exit 1
fi

bcftools view -Ob -v indels "${BAM_BASENAME}.freebayes.bcf" -o "${BAM_BASENAME}.indels_freebayes.bcf"
if [ $? -ne 0 ]; then
    log_message "ERROR" "bcftools view for indels failed"
    exit 1
fi

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

log_message "INFO" "Process completed - SNPs and indels extracted (${DURATION}s)"
log_message "INFO" "Output files: ${BAM_BASENAME}.snps_freebayes.vcf.gz, ${BAM_BASENAME}.indels_freebayes.vcf.gz"

