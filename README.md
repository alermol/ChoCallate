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
- `reference_index_dir`: path to directory with index files for reference genome

After configuration is complete you can run the ChoCallate

```bash
nextflow run main.nf -params-file my_run.yaml
```

The template contains every parameter with its default and a comment that
documents its accepted values, so a run only needs to override the values that
differ. There is no separate command-line help.

Parameters are checked against `nextflow_schema.json` before the first task:
invalid enum values, negative numbers, missing input files and unknown keys
abort with the offending key and the accepted values.

## Resources

Processes use the `process_low`, `process_medium` and `process_high` labels.
Their `cpus`, `memory_mb` and `time` limits are configured in the `resources`
section of the params file. The `mapping.cpu`, `calling.cpu` and
`consensus.cpu` options remain authoritative for their stages and override the
label CPU defaults. JVM heaps are derived from the memory requested for the
task, including GATK processes started by GNU parallel. A task that exits with
status 137 or 140 (typically an out-of-memory kill) is retried up to
`resources.max_retries` times with memory multiplied by the attempt number.
Ordinary caller failures use the same limit.



## Inputs

- **Reads**: FASTQs (`format: "fastq"`) or a pre-aligned BAM (`format: "bam"`). If you provide a BAM, mapping is skipped.
- **Reference genome**: Plain or bgzipped FASTA file. BGZF is detected from the file header; `ref_genome.bgzip` was removed in 4.0.0. Contigs longer than 2,147,483,647 bp (2^31-1) are rejected from the generated `.fai` before mapping starts, because BAM/BCF positions and GATK use 32-bit coordinates.
- **Indexes of reference genome**: All indexes required by the selected mapper/callers. By default the index prefix is the file name of `reference_genome`, including a `.gz` suffix; set `input.reference_index_prefix` when the indexes use a different prefix (for example, when the reference is a symlink). The matching index files are staged as task inputs, so they work with container, cloud and `scratch` execution.
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

Pipeline internals are always BCF. `GENERATE_CONSENSUS` writes one internal
per-sample BCF, `ANNOTATE_HEADER` attaches the command-line header for
per-sample output, `MERGE_OUTPUTS` writes an internal multi-sample BCF when
`output.type` is `"single"`, and `EXPORT_OUTPUT` is the only stage that applies
`output.format`, `remove_invariant` and `split_multiallelic` and publishes the
final file. A VCF output is bgzip-compressed (`.vcf.gz`).

`output.publish_mode` controls how those files are published: `copy` (default)
writes independent copies, while `link` creates symbolic links into the
Nextflow work directory. Linked files require the work directory to be retained
(for example, disable `boost.cleanup`). Results are never moved out of the work
directory, so `-resume` keeps working.


Consensus VCF/BCF records include `FORMAT/GT`, `AD`, and `DP` (AD/DP from BAM), plus site-level `INFO` tags (`AN`, `AC`, `AF`, `NS`, total `DP`, etc.). Details: [Output Structure](https://github.com/alermol/ChoCallate/wiki/Output-Structure).

Consensus rule: every position of the coverage BED is written, monomorphic sites included. Caller records are normalised before voting: multi-allelic records are split into one record per ALT and complex alleles are atomised with `bcftools norm -f ref -m -any --atomize`; each record contributes copies for its own `(REF, ALT)` alleles, so atomised records with different REF alleles at the same POS (for example a SNP `A>G` and a deletion `AT>A`) remain distinct and are all voted on. Split records whose ALT copies exceed the caller's ploidy are inconsistent (for example `1/1` plus `0/1` in a diploid) and make that caller abstain. Alleles are keyed by their REF and ALT sequences, so deletions that share an anchor ALT but have different REF lengths are distinct candidates. Voting has two stages. First, an ALT allele is a candidate when at least `consensus.threshold` callers support it. Second, among the supporters of the winning ALT allele the majority zygosity (heterozygous `0/1` vs homozygous ALT `1/1`) decides the consensus genotype. Ties are resolved by `consensus.tie_policy`: `nocall` (default, no record and independent of caller order) or `first` (the earliest supporter in configured caller order). When several ALT alleles reach the threshold, each becomes its own split consensus record; split records are re-joined for the final output unless `output.split_multiallelic: true`. If no ALT reaches the threshold, the position is reference when at least `consensus.threshold` callers vote reference, otherwise it is a no-call. Callers without a record vote reference by default (`consensus.absent_vote: "ref"`); with `"abstain"` they vote reference only when the sample pileup has at least `coverage.min_coverage` reads at the position. Hom-reference calls (any ploidy) vote reference; caller records with a missing genotype (`./.`, `0/.`) abstain. Caller records whose `FILTER` is not `PASS`/`.` abstain as well by default (`consensus.filtered_vote: "abstain"`); with `"ref"` they vote reference instead. Low-quality records are soft-filtered (`FILTER=LowQual`), not removed, so a caller that saw a variant cannot be mistaken for a reference vote. Positions covered by a deletion are evaluated individually: after the deletion record is written the walk advances by one base, so caller records inside the deletion span are still voted on and positions without a record become reference (`ALT="<*>", GT 0/0`) records when the reference vote reaches the threshold. Overlapping records under the deletion are kept rather than suppressing variant calls inside the span, so `POS..END` covers the deletion and every following BED position is also written.

Reference positions are always emitted internally as gVCF-style records with the symbolic ALT `<*>` (GT `0/0`, `INFO/END`, and `AD=(reference depth, other-allele depth)`). This representation is what lets `bcftools merge --gvcf` map a hom-reference sample's non-reference depth onto another sample's real ALT, so merged `FORMAT/AD` is assigned correctly independently of the requested output format. `output.remove_invariant: true` drops those records only during the final export (per-sample files when `output.type: "sample"`, the merged file otherwise); with the default `false` all covered positions are kept. The final BCF/VCF keeps this `<*>` representation; `output.compress_reference_blocks` only controls whether runs of reference records are collapsed into multi-base blocks.

`output.compress_reference_blocks: true` (default `false`) groups consecutive hom-reference positions without caller records into multi-base gVCF blocks with `INFO/END` during consensus generation, and the blocks are carried through the internal merge into the final output. A block is split at caller records, where depth crosses `coverage.min_coverage`, and wherever BED coverage is not contiguous. `FORMAT/DP` is the minimum depth of the block and `FORMAT/AD` is omitted on multi-base blocks, because one AD value cannot describe every position; when `bcftools merge --gvcf` fills a hom-reference sample at another sample's variant inside such a block, the ALT depth is reported as missing rather than as a fabricated zero. The all-sites representation is preserved because `POS..END` covers every position of the run. Single-base reference records retain their per-position AD, and the option works for both BCF and VCF output.

Calling on fragmented coverage: GATK uses `calling.gatk.interval_padding` (default 150) so indels near the edge of a covered island are assembled, but each padded chunk's calls are restricted to its own unpadded BED interval with `bcftools view -T` before concatenation. Boundary sites therefore come from exactly one chunk, even when the padded neighbours emit different allele representations that `bcftools norm -d exact` cannot deduplicate. Each caller BCF is sorted and deduplicated with `bcftools norm -d exact`. The covered BED is produced by `mosdepth --no-per-base --quantize` (keeping the runs at or above `coverage.min_coverage`, default 5), intersected with the sorted `include_bed`, subtracted with the sorted `exclude_bed`, and split once per sample (`SCATTER_INTERVALS`); the same chunks are reused by every caller and by consensus. Chunk files are concatenated in version order (`ls -v`) and are not re-sorted, except for GATK, which keeps one sort as a guard for padded boundary order. When no base of a sample reaches `coverage.min_coverage`, `GENERATE_COVERAGE` aborts with the sample id and threshold by default. Set `coverage.on_empty: "skip"` to drop that sample with a warning instead; it is then listed in `dropped_samples.txt` and in the merged `SamplesDropped=` header field.

Mapping and BAM input share `bin/postprocess_bam.sh`: reads are coordinate-sorted, MAPQ-filtered after mate fixing and duplicate marking (`bam_filter.min_map_qual`, default 10), optionally deduplicated with `samtools markdup`, and indexed once. Read groups include `PL` (`mapping.platform`, default `ILLUMINA`) and `LB` (`mapping.library`, default the sample id).

Each consensus file carries one structured header line describing the applied rule, for example `##ChoCallate=<Version=4.0.0,Threshold=2,Callers="bcftools,gatk",Rule=two-stage,Tie=nocall,AbsentVote=ref,FilteredVote=abstain,ReferenceBlocks=False>`. `CallersFailed="..."` is added in degraded mode, `SamplesDropped="..."` lists samples that were removed entirely, and `MinCoverage=...` is added when `absent_vote` is `abstain`. The final output also records `##ChoCallateCommand="<full nextflow command>"`. The command line is read from the `NXF_CLI` runtime variable, not embedded in a task script, so changing launch options such as adding `-resume` does not invalidate cached consensus or merge tasks.

Caller selection: `calling.callers` accepts any non-empty subset of `bcftools,freebayes,gatk`; names are trimmed and deduplicated. The configured order only affects the result when `consensus.tie_policy: "first"`; the default `nocall` policy is order-independent. Every caller is retried up to `resources.max_retries` times and a failure then terminates the run. Set `calling.allow_caller_failure: true` to continue without a failed caller instead: a warning is logged, the header records the failed callers, and the run aborts if `consensus.threshold` exceeds the number of surviving callers. If every caller fails for a sample, that sample is dropped with a warning in degraded mode (and with an error otherwise); it is listed in `dropped_samples.txt` and in `SamplesDropped=` in the merged header.

Quality filtering is per caller and per variant type: `calling.<caller>.min_snp_qual` and `calling.<caller>.min_indel_qual` (20/20 in `assets/templates/config.yaml`). Records below the threshold are kept with `FILTER=LowQual` and enter the consensus vote according to `consensus.filtered_vote` (default `abstain`). Both keys are required; the shared `calling.min_snp_qual` key and the `calling.<caller>.min_qual` alias were removed in 4.0.0.

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