#!/bin/bash
set -euo pipefail

export NFX_OFFLINE='true'

echo -e "sample1\t$(readlink -f test_reads_R1.fq.gz)\t$(readlink -f test_reads_R2.fq.gz)" > samples.tsv

cat << EOF > test_config.yaml
input:
  samples_tsv: "$PWD/samples.tsv"
  reference_genome: "$PWD/arth_chr1.fasta.gz"
  reference_index_dir: "$PWD"
  include_bed: null
  exclude_bed: null
  format: "fastq"
  reads_type: "pe"

output:
  directory: "$PWD/ChoCallate_output"
  type: "single"
  format: "bcf"
  remove_invariant: true
  split_multiallelic: false

ref_genome:
  bgzip: true

mapping:
  mapper: "bowtie2"
  extra_args: ""
  cpu: 10

bam_filter:
  min_map_qual: 0

rmdup:
  enabled: false

left_align_indels:
  enabled: false

coverage:
  min_coverage: 2

calling:
  callers: "bcftools,gatk,freebayes"
  min_snp_qual: 20
  cpu: 4
  bcftools:
    call:
      extra_args: ""
    mpileup:
      extra_args: ""
  freebayes:
    extra_args: ""
  gatk:
    extra_args: ""

consensus:
  threshold: 2
  cpu: 4
EOF

echo "Preparing genome index..."
if [[ -e arth_chr1.fasta.gz.1.bt2 ]]; then
  echo "Bowtie2 index already exists, skipping index building."
else
  bowtie2-build --quiet --threads 4 arth_chr1.fasta.gz arth_chr1.fasta.gz
  echo "Bowtie2 index built."
fi
echo "Done"

nextflow run ../main.nf -params-file test_config.yaml

