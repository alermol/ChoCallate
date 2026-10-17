# ChoCallate 🍫

![GitHub Release](https://img.shields.io/github/v/release/alermol/ChoCallate) ![Static Badge](https://img.shields.io/badge/Changelog-orange) ![GitHub License](https://img.shields.io/github/license/alermol/chocallate) ![Static Badge](https://img.shields.io/badge/Wiki-red?link=https%3A%2F%2Fgithub.com%2Falermol%2FChoCallate%2Fwiki)

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
runs the Nextflow entry point to check schema, samplesheet, `--help` and
open-file-limit validation; it needs Nextflow and the pipeline plugins (they
are skipped when `nextflow` is not on `PATH`). Set `CHOCALLATE_NO_ACCEL=1` when
running `generate_consensus.py` directly to force the pure-pysam
implementation.



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

The template contains every parameter with its default, so a run only needs to
override the values that differ. Show the built-in help without any inputs:

```bash
nextflow run main.nf --help            # parameter groups
nextflow run main.nf --help input      # one group in detail
nextflow run main.nf --helpFull        # all parameters
```

Parameters are checked against `nextflow_schema.json` before the first task:
invalid enum values, negative numbers, missing input files and unknown keys
abort with the offending key and the accepted values.



## Inputs

- **Reads**: FASTQs (`format: "fastq"`) or a pre-aligned BAM (`format: "bam"`). If you provide a BAM, mapping is skipped.
- **Reference genome**: Plain or bgzipped FASTA file. BGZF is detected from the file header; `ref_genome.bgzip` was removed in 4.0.0.
- **Indexes of reference genome**: All indexes required by the selected mapper/callers. By default the index prefix is `basename(realpath(reference_genome))`, including a `.gz` suffix; set `input.reference_index_prefix` when the indexes use a different prefix (for example, when the reference is a symlink).
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
- **Tool versions**: `<output.directory>/versions.yml`, written once per run
  from the shared environment (pipeline, mapper, caller, htslib, mosdepth,
  bedops, parallel and Python/pysam versions)

Consensus VCF/BCF records include `FORMAT/GT`, `AD`, and `DP` (AD/DP from BAM), plus site-level `INFO` tags (`AN`, `AC`, `AF`, `NS`, total `DP`, etc.). Details: [Output Structure](https://github.com/alermol/ChoCallate/wiki/Output-Structure).

Consensus rule: every position of the coverage BED is written, monomorphic sites included. Caller records are normalised before voting: multi-allelic records are split into one record per ALT and complex alleles are atomised with `bcftools norm -f ref -m -any --atomize`; all records of one caller at a position form one allele-set vote with remapped GT indices. Voting has two stages. First, an ALT allele is a candidate when at least `consensus.threshold` callers support it. Second, among the supporters of the winning ALT allele the majority zygosity (heterozygous `0/1` vs homozygous ALT `1/1`) decides the consensus genotype. Ties are resolved by `consensus.tie_policy`: `nocall` (default, no record and independent of caller order), `first` (the earliest supporter in configured caller order) or `qual` (largest summed QUAL; a remaining tie is a no-call). When several ALT alleles reach the threshold, each becomes its own split consensus record; split records are re-joined for the final output unless `output.split_multiallelic: true`. If no ALT reaches the threshold, the position is reference when at least `consensus.threshold` callers vote reference, otherwise it is a no-call. Callers without a record vote reference by default (`consensus.absent_vote: "ref"`); with `"abstain"` they vote reference only when the sample pileup has at least `coverage.min_coverage` reads at the position. Hom-reference calls (any ploidy) vote reference; caller records with a missing genotype (`./.`, `0/.`) abstain. Reference positions are emitted as gVCF-style records with the symbolic ALT `<*>` (GT `0/0`, `INFO/END`, and `AD=(reference depth, other-allele depth)`), and merged output uses `bcftools merge --gvcf`, which maps the `<*>` depth onto the real ALT at variant positions of other samples. `output.remove_invariant: true` drops the monomorphic records from the final output (per-sample files when `output.type: "sample"`, the merged file otherwise); with the default `false` all covered positions are kept.

Calling on fragmented coverage: GATK uses `calling.gatk.interval_padding` (default 150) so indels near the edge of a covered island are assembled, and each caller BCF is sorted and deduplicated with `bcftools norm -d exact` to remove records duplicated at chunk boundaries. The covered BED is produced by `mosdepth --no-per-base --quantize` (keeping the runs at or above `coverage.min_coverage`, default 5), intersected with the sorted `include_bed`, subtracted with the sorted `exclude_bed`, and split once per sample (`SCATTER_INTERVALS`); the same chunks are reused by every caller and by consensus. Chunk files are concatenated in version order (`ls -v`) and are not re-sorted, except for GATK, whose interval padding can emit overlapping boundary records.

Mapping and BAM input share `bin/postprocess_bam.sh`: reads are coordinate-sorted, MAPQ-filtered after mate fixing and duplicate marking (`bam_filter.min_map_qual`, default 10), optionally deduplicated with `samtools markdup`, and indexed once. Read groups include `PL` (`mapping.platform`, default `ILLUMINA`) and `LB` (`mapping.library`, default the sample id).

Each consensus file carries one structured header line describing the applied rule, for example `##ChoCallate=<Version=4.0.0,Threshold=2,Callers="bcftools,gatk",Rule=two-stage,Tie=nocall,AbsentVote=ref,RemoveInvariant=False,SplitMultiallelic=False>`. `CallersFailed="..."` is added in degraded mode, and `MinCoverage=...` when `absent_vote` is `abstain`. The same header also records `##ChoCallateCommand="<full nextflow command>"` and `##ChoCallateCallerVersions="bcftools=1.24,..."`; the multi-sample merge re-attaches all `##ChoCallate` lines.

Caller selection: `calling.callers` accepts any non-empty subset of `bcftools,freebayes,gatk`; names are trimmed and deduplicated. The configured order only affects the result when `consensus.tie_policy: "first"`; the default `nocall` policy is order-independent. Every caller is retried twice and a failure then terminates the run. Set `calling.allow_caller_failure: true` to continue without a failed caller instead: a warning is logged, the header records the failed callers, and the run aborts if `consensus.threshold` exceeds the number of surviving callers.

Quality filtering is per caller and per variant type: `calling.<caller>.min_snp_qual` and `calling.<caller>.min_indel_qual` (20/20 in `assets/templates/config.yaml`). Both keys are required; the shared `calling.min_snp_qual` key and the `calling.<caller>.min_qual` alias were removed in 4.0.0.

When `output.type` is `"single"` and `N samples ≥ 2`, single consensus BCF is merged using chunked merge through named pipes. No intermediate BCF files are written to disk. Details: [Merge Strategies and Testing](https://github.com/alermol/ChoCallate/wiki/Merge-Strategies-and-Testing).

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