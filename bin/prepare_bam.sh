#!/usr/bin/env bash
set -euo pipefail

# Args
reads_type="$1"
sample_id="$2"
read1="${3:-}"
read2="${4:-}"
read3="${5:-}"
ref_genome="$6"
genome_dictionary="$7"
genome_fai="$8"
ref_index="$9"
min_map_qual="${10}"
cpus="${11}"
cleanup_intermediate_subfolders="${12}"
cleanup_input_symlinks="${13}"

echo "[$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Process started - Mapping reads (${reads_type})"

echo "[$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Creating intermediate subfolder for BAM preparation"
mkdir -p "${sample_id}_bam_prep"

echo "[$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Running Bowtie2 alignment with ${cpus} threads"

case "${reads_type}" in
  pe)
    bowtie2 --threads ${cpus} --rg-id "${sample_id}" --rg "SM:${sample_id}" -x "${ref_index}" -1 "${read1}" -2 "${read2}" | \
      samtools view -@ ${cpus} -S -b -q ${min_map_qual} -F 4 - | \
      samtools fixmate -@ ${cpus} -m - - | \
      samtools sort -@ ${cpus} -o "${sample_id}_bam_prep/${sample_id}.primary.bam"
    ;;
  se)
    bowtie2 --threads ${cpus} --rg-id "${sample_id}" --rg "SM:${sample_id}" -x "${ref_index}" -U "${read1}" | \
      samtools view -@ ${cpus} -S -b -q ${min_map_qual} -F 4 - | \
      samtools fixmate -@ ${cpus} -m - - | \
      samtools sort -@ ${cpus} -o "${sample_id}_bam_prep/${sample_id}.primary.bam"
    ;;
  mx)
    bowtie2 --threads ${cpus} --rg-id "${sample_id}" --rg "SM:${sample_id}" -x "${ref_index}" -1 "${read1}" -2 "${read2}" -U "${read3}" | \
      samtools view -@ ${cpus} -S -b -q ${min_map_qual} -F 4 - | \
      samtools fixmate -@ ${cpus} -m - - | \
      samtools sort -@ ${cpus} -o "${sample_id}_bam_prep/${sample_id}.primary.bam"
    ;;
  *)
    echo "[$(date -Iseconds)] [ERROR] [BOWTIE2_MAPPING] [${sample_id}] Invalid reads type: ${reads_type}. Available types: se, pe, mx" >&2
    exit 1
    ;;
esac

echo "[$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Running LeftAlignIndels on primary BAM"
gatk LeftAlignIndels -I "${sample_id}_bam_prep/${sample_id}.primary.bam" -O "${sample_id}.bam" -R "${ref_genome}" -OBI false

echo "[$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Indexing final BAM with samtools"
samtools index --csi --threads ${cpus} "${sample_id}.bam"

if [ "${cleanup_intermediate_subfolders}" = "true" ]; then
  rm -rf "${sample_id}_bam_prep"
  echo "[$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Cleaned up intermediate subfolder"
fi

if [ "${cleanup_input_symlinks}" = "true" ]; then
  rm -f "${ref_genome}" "${genome_dictionary}" "${genome_fai}" 2>/dev/null || true
  echo "[$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Cleaned up input files"
fi

echo "[$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Process completed - BAM file created"


