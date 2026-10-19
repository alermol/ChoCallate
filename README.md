# ChoCallate 🍫

[![GitHub Release](https://img.shields.io/github/v/release/alermol/ChoCallate)](https://github.com/alermol/ChoCallate/releases) [![Changelog](https://img.shields.io/badge/Changelog-orange)](CHANGELOG.md) [![GitHub License](https://img.shields.io/github/license/alermol/chocallate)](LICENSE) [![Wiki](https://img.shields.io/badge/Wiki-red)](https://github.com/alermol/ChoCallate/wiki)

**ChoCallate** (**Cho**rus of **Call**ers) - a **Nextflow** pipeline for **consensus-based variant calling**. 

ChoCallate runs several variant callers and applies configurable consensus rules to produce high-confidence **SNVs** and **INDELs**. It addresses a critical challenge in variant calling: individual variant callers can produce different results for the same genomic data, leading to uncertainty in variant identification. By implementing a consensus-driven approach, ChoCallate combines results from multiple state-of-the-art variant callers and applies configurable consensus rules to generate reliable, high-quality variant calls.

## Requirements

- **Linux** (tested). macOS/Windows are not currently tested.
- **Conda** (Miniconda/Anaconda) or **Mamba**
- **Git**
- **Nextflow**

Nextflow resolves the `nf-boost` and `nf-schema` plugins on first use. For an
offline run, install them once on a connected machine with
`nextflow plugin install nf-schema` (and `nf-boost`).



## Installation

```bash
git clone --depth 1 https://github.com/alermol/ChoCallate.git
cd ChoCallate
conda env create -y -f environment.yaml --no-default-packages
conda activate ChoCallate
```

Optionally build the Cython accelerator of `bin/generate_consensus.py` (the
Docker image does this automatically; without it the consensus step falls back
to an equivalent pure-pysam implementation):

```bash
python bin/build_consensus_accel.py
```

Optional verification and cleanup:

```bash
cd test_run
bash run_test.sh
bash cleanup.sh
```

## Testing

The regression suite builds tiny pysam fixtures in code and exercises both the
pure-pysam and the Cython consensus paths, the shared BAM post-processing helper
and the production merge script:

```bash
python bin/build_consensus_accel.py   # optional: enables the Cython-path tests
pytest -q tests/python
```

The suite needs `pytest` (installed by `environment.yaml`) and the environment's
`bcftools`/`samtools`/`pysam`. `tests/python/test_validation.py` additionally
runs the Nextflow entry point to check schema, samplesheet and open-file-limit
validation; it needs Nextflow and the pipeline plugins (they are skipped when
`nextflow` is not on `PATH`). Set `CHOCALLATE_NO_ACCEL=1` when
running `generate_consensus.py` directly to force the pure-pysam
implementation.

The nf-test suite runs the pipeline end to end on the bundled `test_run` data
and checks the all-sites consensus output:

```bash
nf-test test
```



## Docker

ChoCallate is available as a Docker image on DockerHub. For a fuller walkthrough, see the Wiki: [Installing ChoCallate](https://github.com/alermol/ChoCallate/wiki/Installing-ChoCallate#use-the-docker-container).

```bash
docker pull alermol/chocallate:latest
```

Mount your run directory to `/workspace` and run:

```bash
docker run --rm \
  -u "$(id -u):$(id -g)" \
  -v "${PWD}/input_data:/workspace" \
  -w /workspace \
  alermol/chocallate:latest \
  -params-file config.yaml
```

Outputs will be written to the configured `output.directory` (default: `ChoCallate_output`) inside `input_data`.

## Usage

ChoCallate is configured via a Nextflow params YAML file. Start from the template.

```bash
cp assets/templates/config.yaml my_run.yaml
```

Minimum set of parameters in `my_run.yaml`:

- `samples_tsv`: input samples TSV (formats below)
- `reference_genome`: reference FASTA (plain or bgzip-compressed)
- `reference_index_dir`: path to directory with index files for reference genome (required for `format: fastq`, unused for `format: bam`)

After configuration is complete you can run the ChoCallate

```bash
nextflow run main.nf -params-file my_run.yaml
```

The template contains every parameter with its default and a comment that
documents its accepted values, so a run only needs to override the values that
differ. There is no separate command-line help.

Parameters are checked against `nextflow_schema.json` before the first task:
invalid enum values, values below the accepted minimum, missing input files
and unknown keys abort with the offending key and the accepted values.

`boost.cleanup` (default `"v2"`) controls the nf-boost work-directory
cleanup. The cleanup observers delete the files of a completed task once no
downstream process needs them, and deleted files cannot be reused. A later
`nextflow run ... -resume` therefore re-executes every task whose files were
removed. Set `boost.cleanup: false` in the params file to keep the work
directory and let `-resume` reuse completed tasks; `"v1"` and `"v2"` enable
the corresponding cleanup observer.

## Resources

Processes use the `process_low`, `process_medium` and `process_high` labels.
Their `cpus`, `memory_mb` and `time` limits are configured in the `resources`
section of the params file. The `mapping.cpu`, `calling.cpu` and
`consensus.cpu` options remain authoritative for their stages and override the
label CPU defaults. JVM heaps are derived from the memory requested for the
task, including GATK processes started by GNU parallel. A task that exits with
status 137 (typically an out-of-memory kill) or 140 (a walltime kill) is
retried up to `resources.max_retries` times with memory and time multiplied by
the attempt number. `resources.limits` (`cpus`, `memory_mb`, `time`) caps the
escalated request, so a late retry cannot ask for resources the executor will
never provide. A caller failure that is not transient (for example an
invalid extra argument) terminates immediately; set
`calling.allow_caller_failure: true` to ignore it and continue with the
surviving callers. GNU parallel
is started with `--halt now,fail=1` and its job log is inspected, so a
SIGKILLed child (signal 9, exit 137, or the GATK launcher's 247) is reported
as task exit 137 instead of a generic non-zero parallel status. Mapping tasks
split `mapping.cpu` between the mapper and the concurrent BAM
post-processing. Because `samtools sort` runs in the same task while the
mapper index is resident, it receives half of the task memory that remains
after the aligner reserve, split across its threads. The reserve is
`mapping.aligner_reserve_mb`; when it is 0 (the default), the pipeline
estimates it from the staged index files plus 1 GiB. For large genomes the
index dominates the task memory: leave the estimate in place or set the
reserve explicitly (roughly 4 GiB for a human bowtie2 index, 6 GiB for BWA
and 8 GiB or more for a minimap2 `.mmi`) and raise
`resources.medium.memory_mb` so that the reserve, the sort buffer and the
streaming stages fit together. A reserve that is not smaller than the task
memory stops the task with a message naming both values.



## Inputs

- **Reads**: FASTQs (`format: "fastq"`) or a pre-aligned BAM (`format: "bam"`). If you provide a BAM, mapping is skipped.
- **Reference genome**: Plain or bgzipped FASTA file. BGZF is detected from the file header; `ref_genome.bgzip` was removed in 4.0.0. Contigs longer than 2,147,483,647 bp (2^31-1) are rejected from the generated `.fai` before mapping starts, because BAM/BCF positions and GATK use 32-bit coordinates.
- **Indexes of reference genome**: All indexes required by the selected mapper. By default the index prefix is the file name of `reference_genome`, including a `.gz` suffix; set `input.reference_index_prefix` when the indexes use a different prefix (for example, when the reference is a symlink). The prefix is a plain file name in `reference_index_dir` and cannot contain path separators, because indexes are staged under a flat `index/` directory. A complete mapper index is staged as a task input: the six `.bt2` or six `.bt2l` files for bowtie2, `.amb`/`.ann`/`.bwt`/`.pac`/`.sa` for bwa, or the `.mmi` for minimap2. A minimap2 index must be built with `minimap2 -x sr -d`: the k-mer size and window stored in the `.mmi` override the `-x sr` preset at mapping time, and minimap2 only emits a warning. An incomplete set aborts before the first task and lists the missing files. Staged indexes work with container, cloud and `scratch` execution.
- **BED regions**: `include_bed` and `exclude_bed` may be unsorted and may contain `track`, `browser`, and `#` comment lines; these lines are removed before sorting. Every contig must be present in the reference `.fai`; otherwise `GENERATE_COVERAGE` fails with the BED file and the offending contig names instead of reporting empty coverage.
- **Paths**: Read paths in the samplesheet may be absolute or relative to the samplesheet. Reference, index and BED paths are resolved relative to the launch directory; absolute paths are recommended for those.

`samples_tsv` formats:

- **FASTQ + paired-end** (`reads_type: "pe"`): `sample_id<TAB>R1<TAB>R2`
- **FASTQ + single-end** (`reads_type: "se"`): `sample_id<TAB>R1`
- **FASTQ + mixed** (`reads_type: "mx"`): `sample_id<TAB>R1<TAB>R2<TAB>U` (Bowtie2 mapping only)
- **BAM** (`format: "bam"`): `sample_id<TAB>bam_path`

Blank lines and lines starting with `#` are ignored. An optional header row
(`sample_id`, `sample`, `sampleid` or `id` in the first column) is skipped.
Sample ids must be unique and may contain letters, digits, `.`, `_` and `-`
(must start with a letter or digit). Every data row must have exactly the
number of columns shown above, and every referenced read file must exist.



## Outputs

Published outputs are written to `output.directory` (default: `ChoCallate_output`), including standard Nextflow reports:

- `pipeline_report.html`
- `timeline_report.html`
- `trace.txt`

Consensus outputs depend on `output.type` and `output.format`:

- **Single merged** (default): `<output.directory>/consensus.bcf` or `<output.directory>/consensus.vcf.gz`
- **Per-sample**: `<output.directory>/per_sample/<sample_id>.bcf` or `<sample_id>.vcf.gz`

`dropped_samples.txt` is always published in `output.directory`. It lists,
one per line, the samples that were removed from the run: samples whose
every caller failed as well as samples with an empty coverage BED when
`coverage.on_empty: "skip"` is set.

Per-sample mapping statistics are published under
`<output.directory>/mapping_stats/` for every sample that reached mapping:

- `<sample_id>.raw.flagstat`: `samtools flagstat` of the aligner output
  before MAPQ filtering and duplicate handling. For BAM input this is the
  input BAM. It provides the mapping rate that the post-processing filters
  would otherwise hide.
- `<sample_id>.final.flagstat`, `<sample_id>.stats` and
  `<sample_id>.idxstats`: `samtools flagstat`, `samtools stats` and
  `samtools idxstats` of the final BAM consumed by coverage and calling.
- `mapping_summary.tsv`: one row per samplesheet entry, in samplesheet order,
  with alignment-record counts, raw mapped and properly paired fractions,
  duplicate and MQ0 fractions, mean read length and MAPQ, insert size, error
  rate and the number of contigs with mapped reads.

Counts are alignment records: a paired-end read pair contributes two records.
The raw report is taken before `bam_filter.min_map_qual` and `rmdup.mode` are
applied, so `raw_mapped_pct` describes the mapper output rather than the
filtered BAM. Statistics are reported even for samples that are later dropped
by the coverage filter or a failed caller.

Pipeline internals are always BCF. `GENERATE_CONSENSUS` writes one internal
per-sample BCF, `MERGE_OUTPUTS` writes an internal multi-sample BCF when
`output.type` is `"single"`, and `EXPORT_OUTPUT` is the only stage that applies
`output.format`, `remove_invariant` and `split_multiallelic` and publishes the
final file. A VCF output is bgzip-compressed (`.vcf.gz`).

`output.publish_mode` controls how those files are published: `copy` (default)
writes independent copies, while `link` creates symbolic links into the
Nextflow work directory. Linked files require the work directory to be retained
(for example, with `boost.cleanup: false`). Results are never moved out of the
work directory, so `-resume` keeps working.


Consensus VCF/BCF records include `FORMAT/GT`, `AD`, and `DP` (AD/DP from BAM), plus site-level `INFO` tags (`AN`, `AC`, `AF`, `NS`, total `DP`, etc.). Details: [Output Structure](https://github.com/alermol/ChoCallate/wiki/Output-Structure).

Consensus rule: every position of the coverage BED is evaluated. Monomorphic sites are written as reference records unless `output.remove_invariant: true`; with `output.type: "sample"` and `output.remove_invariant: true`, reference positions are skipped during consensus because the export would drop them. Positions whose reference base is an IUPAC ambiguity code (N included) and positions that end as a consensus no-call are not written: the first cannot be represented with an unambiguous REF allele, and the second has no consensus genotype to report. Caller records are normalised before voting: multi-allelic records are split into one record per ALT and complex alleles are atomised with `bcftools norm -f ref -m -any --atomize`; each record contributes copies for its own `(REF, ALT)` alleles, so atomised records with different REF alleles at the same POS (for example a SNP `A>G` and a deletion `AT>A`) remain distinct and are all voted on. Split records whose ALT copies exceed the caller's ploidy are inconsistent (for example `1/1` plus `0/1` in a diploid) and make that caller abstain. Alleles are keyed by their REF and ALT sequences, so deletions that share an anchor ALT but have different REF lengths are distinct candidates. Voting has two stages. First, an ALT allele is a candidate when at least `consensus.threshold` callers support it. Second, among the supporters of the winning ALT allele the majority zygosity (heterozygous `0/1` vs homozygous ALT `1/1`) decides the consensus genotype. Ties are resolved by `consensus.tie_policy`: `nocall` (default, no record and independent of caller order) or `first` (the earliest supporter in configured caller order). When several ALT alleles reach the threshold, each becomes its own split consensus record; split records are re-joined for the final output unless `output.split_multiallelic: true`. The combined ALT copies of those records must fit the ploidy: if they do not (for example two hom-ALT alleles at a diploid site), the position is a no-call instead of producing split records that cannot be represented by one genotype. If no ALT reaches the threshold, the position is reference when at least `consensus.threshold` callers vote reference, otherwise it is a no-call. Callers without a record vote reference by default (`consensus.absent_vote: "ref"`); with `"abstain"` they vote reference only when the sample pileup has at least `coverage.min_coverage` reads at the position. Hom-reference calls (any ploidy) vote reference; caller records with a missing genotype (`./.`, `0/.`) abstain. Caller records whose `FILTER` is not `PASS`/`.` abstain as well by default (`consensus.filtered_vote: "abstain"`); with `"ref"` they vote reference instead. Low-quality records are soft-filtered (`FILTER=LowQual`), not removed, so a caller that saw a variant cannot be mistaken for a reference vote. Positions covered by a deletion are evaluated individually: after the deletion record is written the walk advances by one base, so caller records inside the deletion span are still voted on and variant calls under the deletion are kept. Reference records are not written over a called deletion's REF span, because a `0/0` call over deleted bases would contradict the deletion genotype; those positions are represented by the deletion record's `POS..END`, and reference output resumes after the span.

Reference positions are emitted internally as gVCF-style records with the symbolic ALT `<*>` (GT `0/0`, `INFO/END`, and `AD=(reference depth, other-allele depth)`), except when `output.type: "sample"` and `output.remove_invariant: true`: per-sample output has no merge step, and the export would drop those records anyway, so the consensus stage is told to write variants only. In the merged path this representation is what lets `bcftools merge --gvcf` map a hom-reference sample's non-reference depth onto another sample's real ALT, so merged `FORMAT/AD` is assigned correctly independently of the requested output format. `output.remove_invariant: true` drops reference records during the final export; with the default `false`, every covered position that has a reference or variant consensus is kept, except the no-call and IUPAC-ambiguous positions described above. The final BCF/VCF keeps this `<*>` representation; `output.compress_reference_blocks` only controls whether runs of reference records are collapsed into multi-base blocks.

`output.compress_reference_blocks: true` is the default: it groups consecutive hom-reference positions without caller records into multi-base gVCF blocks with `INFO/END` during consensus generation, and the blocks are carried through the internal merge into the final output. A block is split at caller records, where depth crosses `coverage.min_coverage`, and wherever BED coverage is not contiguous. `FORMAT/DP` is the minimum depth of the block and `FORMAT/AD` is omitted on multi-base blocks, because one AD value cannot describe every position; when `bcftools merge --gvcf` fills a hom-reference sample at another sample's variant inside such a block, the ALT depth is reported as missing rather than as a fabricated zero. The all-sites representation is preserved because `POS..END` covers every position of the run. Single-base reference records retain their per-position AD, and the option works for both BCF and VCF output.

Calling on fragmented coverage: GATK uses `calling.gatk.interval_padding` (default 150) so indels near the edge of a covered island are assembled, but each padded chunk's calls are restricted to its own unpadded BED interval with `bcftools view -T` before concatenation. Boundary sites therefore come from exactly one chunk, even when the padded neighbours emit different allele representations that `bcftools norm -d exact` cannot deduplicate. Each caller BCF is sorted and deduplicated with `bcftools norm -d exact`. The covered BED is produced by `mosdepth --no-per-base --quantize` (keeping the runs at or above `coverage.min_coverage`, default 5), intersected with the sorted `include_bed`, subtracted with the sorted `exclude_bed`, and split once per sample (`SCATTER_INTERVALS`); the same chunks are reused by every caller and by consensus. Chunk files are concatenated in version order (`ls -v`) and are not re-sorted, except for GATK, which keeps one sort as a guard for padded boundary order. When no base of a sample reaches `coverage.min_coverage`, `GENERATE_COVERAGE` aborts with the sample id and threshold by default. Set `coverage.on_empty: "skip"` to drop that sample with a warning instead; it is then listed in `dropped_samples.txt` and in the merged `SamplesDropped=` header field.

Mapping and BAM input share `bin/postprocess_bam.sh`: aligner output is piped through `samtools view -u` as uncompressed BAM, then reads are coordinate-sorted and MAPQ-filtered (`bam_filter.min_map_qual`, default 10) after mate fixing and duplicate handling; when duplicate handling is disabled, the MAPQ filter runs before sorting. The post-processing CPU allocation is divided among its concurrent samtools stages, and each BAM is indexed once. `rmdup.mode` selects duplicate handling with `samtools markdup`: `none` leaves duplicates unflagged, `mark` (default) flags them without removing them, and `remove` flags and removes them from the BAM. The callers, mosdepth coverage and the consensus pileup already skip duplicate-flagged reads, so `mark` is sufficient for calling while keeping the reads available for inspection. The former `rmdup.enabled` key was replaced by `rmdup.mode`. Read groups include `PL` (`mapping.platform`, default `ILLUMINA`) and `LB` (`mapping.library`, default the sample id). Both values accept only letters, digits, `_`, `.`, `:` and `-`.

Each consensus file carries one structured header line describing the applied rule, for example `##ChoCallate=<Version=4.0.0,Threshold=2,Callers="bcftools,gatk",Rule=two-stage,Tie=nocall,AbsentVote=ref,FilteredVote=abstain,ReferenceBlocks=True>`. `CallersFailed="..."` is added in degraded mode, `SamplesDropped="..."` lists samples that were removed entirely, and `MinCoverage=...` is added when `absent_vote` is `abstain`. In a merged multi-sample file, `Callers` and `CallersFailed` are the union across samples; when the per-sample values differ (for example because different callers failed in degraded mode), each sample's own values are recorded in a `##ChoCallateSamples=<Sample=...,Callers="...",CallersFailed="...">` line.

Caller selection: `calling.callers` accepts any non-empty subset of `bcftools,freebayes,gatk`; names are trimmed and deduplicated. The configured order only affects the result when `consensus.tie_policy: "first"`; the default `nocall` policy is order-independent. Every caller retries transient failures (exit status 137 or 140) up to `resources.max_retries` times; any other failure terminates the run immediately. Set `calling.allow_caller_failure: true` to continue without a failed caller instead: a warning is logged, the header records the failed callers, and the run aborts if `consensus.threshold` exceeds the number of surviving callers. If every caller fails for a sample, that sample is dropped with a warning in degraded mode (and with an error otherwise); it is listed in `dropped_samples.txt` and in `SamplesDropped=` in the merged header.

Quality filtering is per caller and per variant type: `calling.<caller>.min_snp_qual` and `calling.<caller>.min_indel_qual` (20/20 in `assets/templates/config.yaml`). The filter runs after `bcftools norm -m -any --atomize`, so multi-allelic, MNP and complex records are split first and each resulting record is filtered by its own type (an atomised SNP uses `min_snp_qual`, not the indel threshold). Records below the threshold are kept with `FILTER=LowQual` and enter the consensus vote according to `consensus.filtered_vote` (default `abstain`). Both keys are required; the shared `calling.min_snp_qual` key and the `calling.<caller>.min_qual` alias were removed in 4.0.0.

Caller defaults are overridden explicitly so they cannot silently limit calling. `bcftools mpileup` receives `--max-depth` and `--max-idepth` from `calling.bcftools.mpileup.max_depth` and `calling.bcftools.mpileup.max_idepth` (both default 100000) instead of the bcftools defaults of 250, which down-sample high-coverage libraries for SNP and indel calling. GATK HaplotypeCaller receives `--standard-min-confidence-threshold-for-calling min(min_snp_qual, min_indel_qual)`, so thresholds below the GATK default of 30 take effect. FreeBayes receives `--use-best-n-alleles` from `calling.freebayes.use_best_n_alleles` (default 4; 0 evaluates all alleles) and `--limit-coverage` from `calling.freebayes.limit_coverage` (default 100000), which bound allele enumeration and per-sample coverage in repeats and high-depth regions.

When `output.type` is `"single"`, per-sample consensus BCFs are merged into one internal BCF using chunked merge through named pipes. No intermediate merge files are written to disk; the result is converted by `EXPORT_OUTPUT`. Details: [Merge Strategies and Testing](https://github.com/alermol/ChoCallate/wiki/Merge-Strategies-and-Testing).

## Additional documentation

- **Wiki home**: [ChoCallate Wiki](https://github.com/alermol/ChoCallate/wiki)
- **Install**: [Installing ChoCallate](https://github.com/alermol/ChoCallate/wiki/Installing-ChoCallate)
- **Quick start / config & CLI**: [CLI Reference](https://github.com/alermol/ChoCallate/wiki/CLI-Reference)
- **Multi-sample merge**: [Merge Strategies and Testing](https://github.com/alermol/ChoCallate/wiki/Merge-Strategies-and-Testing)
- **Outputs**: [Output Structure](https://github.com/alermol/ChoCallate/wiki/Output-Structure)



## Contribution

See [CONTRIBUTING.md](https://github.com/alermol/ChoCallate/blob/main/CONTRIBUTING.md)

## License

[MIT](https://github.com/alermol/ChoCallate/blob/main/LICENSE)

## Development roadmap

See [Development Roadmap](https://github.com/alermol/ChoCallate/blob/main/ROADMAP.md) for planned container support and additional callers/mappers.