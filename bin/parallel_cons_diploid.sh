#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

rc=$1 # region coordinates
rn=$2 # region number
sample_name=$3
mutation_type=$4


region=$(echo $rc | sed 's/ /:/' | sed 's/ /-/')

for i in $(ls ${sample_name}.${mutation_type}_* | grep -v '.csi')
do
    bcftools view -r ${region} ${i} > ${rn}.${i}
done

bcftools view -r ${region} zero.vcf.gz > ${rn}.zero.vcf

if [ ${mutation_type} == "snps" ]; then
    process_snps.py \
        --zero_vcf ${rn}.zero.vcf \
        --vcf1 ${rn}.${sample_name}.snps_bcftools \
        --vcf2 ${rn}.${sample_name}.snps_freebayes \
        --vcf3 ${rn}.${sample_name}.snps_gatk \
        --vcf4 ${rn}.${sample_name}.snps_vardict \
        --vcf5 ${rn}.${sample_name}.snps_snver \
        --sample ${sample_name} \
        --chr ${rn}
else
    process_indels.py \
        --vcf1 ${rn}.${sample_name}.indels_bcftools \
        --vcf2 ${rn}.${sample_name}.indels_freebayes \
        --vcf3 ${rn}.${sample_name}.indels_gatk \
        --vcf4 ${rn}.${sample_name}.indels_vardict \
        --vcf5 ${rn}.${sample_name}.indels_snver \
        --sample ${sample_name} \
        --chr ${rn}
fi


for i in $(ls ${sample_name}.${mutation_type}_* | grep -v '.csi')
do
    rm ${rn}.${i}
done

rm ${rn}.zero.vcf

find all_chrs/ -name '*.vcf' -type f -exec bgzip {} \;
