#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

rc=$1 # region coordinates
rn=$2 # region number
sample_name=$3
mutation_type=$4
cons_threshold=$5

region=$(echo $rc | sed 's/ /:/' | sed 's/ /-/')

# Logging function
log_message() {
    local level=$1
    local message=$2
    local timestamp=$(date -Iseconds)
    echo "[${timestamp}] [${level}] [CONSENSUS_GENERATION] [${sample_name}:${rn}] ${message}"
}

# Log process start
log_message "INFO" "Process started - Consensus generation for ${mutation_type}"
log_message "INFO" "Parameters: region=${region}, threshold=${cons_threshold}, sample=${sample_name}"

# Record start time
START_TIME=$(date +%s)

# Extract region-specific VCFs
log_message "INFO" "Extracting region-specific VCFs"
for i in $(ls ${sample_name}.${mutation_type}_* | grep -v '.csi')
do
    bcftools view -r ${region} ${i} > ${rn}.${i}
    log_message "DEBUG" "Extracted region ${region} from ${i}"
done

# Extract region-specific zero VCF
log_message "INFO" "Extracting region-specific zero VCF"
bcftools view -r ${region} zero.vcf.gz > ${rn}.zero.vcf

# Process based on mutation type
if [ ${mutation_type} == "snps" ]; then
    log_message "INFO" "Processing SNPs consensus"
    process_snps.py \
        --zero_vcf ${rn}.zero.vcf \
        --vcfs $(ls ${rn}.${sample_name}.snps_*.vcf.gz | tr '\n' ',') \
        --sample ${sample_name} \
        --chr ${rn} \
        --cons_threshold ${cons_threshold}
    
    if [ $? -eq 0 ]; then
        log_message "INFO" "SNPs consensus processing completed successfully"
    else
        log_message "ERROR" "SNPs consensus processing failed"
        exit 1
    fi
else
    log_message "INFO" "Processing indels consensus"
    process_indels.py \
        --vcfs $(ls ${rn}.${sample_name}.indels_*.vcf.gz | tr '\n' ',') \
        --sample ${sample_name} \
        --chr ${rn} \
        --cons_threshold ${cons_threshold}
    
    if [ $? -eq 0 ]; then
        log_message "INFO" "Indels consensus processing completed successfully"
    else
        log_message "ERROR" "Indels consensus processing failed"
        exit 1
    fi
fi

# Cleanup temporary files
log_message "INFO" "Cleaning up temporary files"
for i in $(ls ${sample_name}.${mutation_type}_* | grep -v '.csi')
do
    rm ${rn}.${i}
done

rm ${rn}.zero.vcf

# Compress output VCFs
log_message "INFO" "Compressing output VCFs"
find all_chrs/ -name '*.vcf' -type f -exec bgzip {} \;

# Calculate duration and log completion
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

log_message "INFO" "Process completed - Consensus generated for region ${region}"
log_message "INFO" "Performance: ${DURATION} seconds"
