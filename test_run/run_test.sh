#!/bin/bash
set -euo pipefail

export NXF_OFFLINE='true'

echo -e "sample1\t$(readlink -f test_reads_R1.fq.gz)\t$(readlink -f test_reads_R2.fq.gz)" > samples.tsv

cat << EOF > test_config.yaml
input:
  samples_tsv: "$PWD/samples.tsv"
  reference_genome: "$PWD/arth_chr1.fasta.gz"
  reference_index_dir: "$PWD"
  reference_index_prefix: null
  include_bed: null
  exclude_bed: null
  format: "fastq"
  reads_type: "pe"

output:
  directory: "$PWD/ChoCallate_output"
  type: "single"
  format: "bcf"
  remove_invariant: false
  split_multiallelic: false

mapping:
  mapper: "bowtie2"
  extra_args: ""
  cpu: 10
  platform: "ILLUMINA"
  library: null

bam_filter:
  min_map_qual: 10

rmdup:
  mode: "none"

left_align_indels:
  enabled: false

coverage:
  min_coverage: 5

calling:
  callers: "bcftools,gatk,freebayes"
  cpu: 4
  allow_caller_failure: false
  bcftools:
    min_snp_qual: 20
    min_indel_qual: 20
    call:
      extra_args: ""
    mpileup:
      extra_args: ""
  freebayes:
    min_snp_qual: 20
    min_indel_qual: 20
    extra_args: ""
  gatk:
    min_snp_qual: 20
    min_indel_qual: 20
    interval_padding: 150
    extra_args: ""

consensus:
  threshold: 2
  cpu: 4
  absent_vote: "ref"
  tie_policy: "nocall"
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

