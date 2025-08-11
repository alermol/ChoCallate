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

process_${mutation_type}.py --vcf1 ${rn}.${sample_name}.${mutation_type}_bcftools --vcf2 ${rn}.${sample_name}.${mutation_type}_freebayes --vcf3 ${rn}.${sample_name}.${mutation_type}_gatk --vcf4 ${rn}.${sample_name}.${mutation_type}_vardict --vcf5 ${rn}.${sample_name}.${mutation_type}_snver --sample ${sample_name} --chr ${rn}


for i in $(ls ${sample_name}.${mutation_type}_* | grep -v '.csi')
do
    rm ${rn}.${i}
done

find all_chrs/ -name '*.vcf' -type f -exec bgzip {} \;
