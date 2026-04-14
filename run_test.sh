#!/bin/bash
set -euo pipefail

echo "Preparing genome index..."
if [[ -e test_data/arth_chr1.fasta.gz.1.bt2 ]]; then
  echo "Bowtie2 index already exists, skipping index building."
else
  bowtie2-build --quiet --threads 4 test_data/arth_chr1.fasta.gz test_data/arth_chr1.fasta.gz
  echo "Bowtie2 index built."
fi
echo "Done"

name="sample1"
read1=$(readlink -f test_data/test_reads_R1.fq.gz)
read2=$(readlink -f test_data/test_reads_R2.fq.gz)
echo -e "sample1\t${read1}\t${read2}" > test_data/samples.tsv

nextflow run main.nf -params-file test_data/test_config.yaml

