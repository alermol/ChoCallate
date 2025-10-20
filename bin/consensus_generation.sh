#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

rc=$1 # region coordinates
rn=$2 # region number
sample_name=$3
mutation_type=$4
cons_threshold=$5

region=$(echo $rc | sed 's/ /:/' | sed 's/ /-/')

log_message() {
    local level=$1
    local message=$2
    local timestamp=$(date -Iseconds)
    echo "[${timestamp}] [${level}] [CONSENSUS_GENERATION] [${sample_name}:${rn}] ${message}"
}

log_message "INFO" "Process started - Consensus generation for ${mutation_type} (threshold=${cons_threshold})"

START_TIME=$(date +%s)

log_message "INFO" "Extracting region-specific BCFs for region ${region}"
for i in $(ls ${sample_name}.${mutation_type}_* | grep -v '.csi')
do
    bcftools view -r ${region} ${i} > ${rn}.${i}
done

log_message "INFO" "Extracting region-specific zero BCF"
bcftools view -r ${region} zero.bcf > ${rn}.zero.bcf

log_message "INFO" "Extracting chromosome names from zero BCF"
chr_names=$(cat zero_chr_names.txt)

if [ ${mutation_type} == "snps" ]; then
    log_message "INFO" "Processing SNPs consensus"
    process_snps.py \
        --zero_bcf ${rn}.zero.bcf \
        --bcfs "$(ls ${rn}.${sample_name}.snps_*.bcf | tr '\n' ',')" \
        --sample ${sample_name} \
        --chr "${chr_names}" \
        --cons_threshold ${cons_threshold} | bcftools view -Ob -o all_chrs/${rn}.bcf

    if [ $? -eq 0 ]; then
        log_message "INFO" "SNPs consensus processing completed successfully"
    else
        log_message "ERROR" "SNPs consensus processing failed"
        exit 1
    fi
else
    log_message "INFO" "Processing indels consensus"
    process_indels.py \
        --zero_bcf ${rn}.zero.bcf \
        --bcfs "$(ls ${rn}.${sample_name}.indels_*.bcf | tr '\n' ',')" \
        --sample ${sample_name} \
        --chr "${chr_names}" \
        --cons_threshold ${cons_threshold} | bcftools view -Ob -o all_chrs/${rn}.bcf
    
    if [ $? -eq 0 ]; then
        log_message "INFO" "Indels consensus processing completed successfully"
    else
        log_message "ERROR" "Indels consensus processing failed"
        exit 1
    fi
fi

log_message "INFO" "Cleaning up temporary files"
for i in $(ls ${sample_name}.${mutation_type}_* | grep -v '.csi')
do
    if [ -e "${rn}.${i}" ]; then
        rm "${rn}.${i}"
    fi
done

if [ -e "${rn}.zero.bcf" ]; then
    rm -f "${rn}.zero.bcf"
fi

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

log_message "INFO" "Process completed - Consensus generated for region ${region} (${DURATION}s)"
