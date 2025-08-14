#!/bin/bash

# FreeBayes variant caller script

set -e

# Input parameters
BAM_FILE=$1
COVERAGE_FILE=$2
REF_GENOME=$3
PLOIDY=$4
READS_SOURCE=$5
MIN_BASE_QUALITY=$6
MIN_SNP_QUAL=$7

# Extract base name from BAM file
BAM_BASENAME=$(basename "$BAM_FILE" .bam)

# Logging function
log_message() {
    local level=$1
    local message=$2
    local timestamp=$(date -Iseconds)
    echo "[${timestamp}] [${level}] [FREEBAYES_CALLER] [${BAM_BASENAME}] ${message}"
}

# Log process start
log_message "INFO" "Process started - FreeBayes variant calling"
log_message "INFO" "Parameters: BAM=${BAM_FILE}, REF=${REF_GENOME}, PLOIDY=${PLOIDY}, READS_SOURCE=${READS_SOURCE}, MIN_BQ=${MIN_BASE_QUALITY}, MIN_QUAL=${MIN_SNP_QUAL}"

# Record start time
START_TIME=$(date +%s)

# Run FreeBayes with appropriate parameters based on reads source
log_message "INFO" "Running FreeBayes variant calling with ${READS_SOURCE} parameters"
if [[ "$READS_SOURCE" == "gbs" ]]; then
    freebayes --fasta-reference "$REF_GENOME" --targets "$COVERAGE_FILE" --dont-left-align-indels --ploidy "$PLOIDY" \
        --use-best-n-alleles 4 --min-alternate-qsum "$MIN_BASE_QUALITY" --hwe-priors-off --no-population-priors \
        --binomial-obs-priors-off --allele-balance-priors-off --min-base-quality "$MIN_BASE_QUALITY" \
        --haplotype-length -1 --throw-away-complex-obs --no-partial-observations --bam "$BAM_FILE" --limit-coverage 250 | \
        bcftools filter -e"QUAL<$MIN_SNP_QUAL" - | \
        bcftools view --min-alleles 2 --max-alleles 2 - | bcftools annotate --force -x INFO,FORMAT - | \
        bcftools norm --fasta-ref "$REF_GENOME" --atom-overlaps '.' --atomize -Oz -o "${BAM_BASENAME}.freebayes.vcf.gz"
else
    freebayes --fasta-reference "$REF_GENOME" --targets "$COVERAGE_FILE" --dont-left-align-indels --ploidy "$PLOIDY" \
        --use-best-n-alleles 4 --min-alternate-qsum "$MIN_BASE_QUALITY" --hwe-priors-off --no-population-priors \
        --allele-balance-priors-off --min-base-quality "$MIN_BASE_QUALITY" \
        --haplotype-length -1 --throw-away-complex-obs --no-partial-observations --bam "$BAM_FILE" --limit-coverage 250 | \
        bcftools filter -e"QUAL<$MIN_SNP_QUAL" - | \
        bcftools view --min-alleles 2 --max-alleles 2 - | bcftools annotate --force -x INFO,FORMAT - | \
        bcftools norm --fasta-ref "$REF_GENOME" --atom-overlaps '.' --atomize -Oz -o "${BAM_BASENAME}.freebayes.vcf.gz"
fi

if [ $? -eq 0 ]; then
    log_message "INFO" "FreeBayes variant calling completed successfully"
else
    log_message "ERROR" "FreeBayes variant calling failed"
    exit 1
fi

# Extract SNPs and indels
log_message "INFO" "Extracting SNPs and indels from FreeBayes output"
bcftools view -v snps -Oz -o "${BAM_BASENAME}.snps_freebayes.vcf.gz" "${BAM_BASENAME}.freebayes.vcf.gz"
bcftools view -v indels -Oz -o "${BAM_BASENAME}.indels_freebayes.vcf.gz" "${BAM_BASENAME}.freebayes.vcf.gz"

# Calculate duration and log completion
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

log_message "INFO" "Process completed - SNPs and indels extracted"
log_message "INFO" "Performance: ${DURATION} seconds"
log_message "INFO" "Output files: ${BAM_BASENAME}.snps_freebayes.vcf.gz, ${BAM_BASENAME}.indels_freebayes.vcf.gz"

