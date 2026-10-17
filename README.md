# ChoCallate 🍫

![GitHub Release](https://img.shields.io/github/v/release/alermol/ChoCallate) ![Static Badge](https://img.shields.io/badge/Changelog-orange) ![GitHub License](https://img.shields.io/github/license/alermol/chocallate) ![Static Badge](https://img.shields.io/badge/Wiki-red?link=https%3A%2F%2Fgithub.com%2Falermol%2FChoCallate%2Fwiki)

**ChoCallate** (**Cho**rus of **Call**ers) - a **Nextflow** pipeline for **consensus-based variant calling**. 

ChoCallate runs several variant callers and applies configurable consensus rules to produce high-confidence **SNVs** and **INDELs**. It addresses a critical challenge in variant calling: individual variant callers can produce different results for the same genomic data, leading to uncertainty in variant identification. By implementing a consensus-driven approach, ChoCallate combines results from multiple state-of-the-art variant callers and applies configurable consensus rules to generate reliable, high-quality variant calls.

## Requirements

- **Linux** (tested). macOS/Windows are not currently tested.
- **Conda** (Miniconda/Anaconda) or **Mamba**
- **Git**
- **Nextflow**



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
pure-pysam and the Cython consensus paths as well as the merge chain:

```bash
python bin/build_consensus_accel.py   # optional: enables the Cython-path tests
pytest -q tests/python
```

The suite needs `pytest` (installed by `environment.yaml`) and the environment's
`bcftools`/`pysam`. Set `CHOCALLATE_NO_ACCEL=1` when running
`generate_consensus.py` directly to force the pure-pysam implementation.



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

Outputs will be written to the configured `outdir` (default: `ChoCallate_output`) inside `input_data`.

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



## Inputs

- **Reads**: FASTQs (`format: "fastq"`) or a pre-aligned BAM (`format: "bam"`). If you provide a BAM, mapping is skipped.
- **Reference genome**: Plain or bgzipped FASTA file.
- **Indexes of reference genome**: All indexes required by the selected mapper/callers.
- **Tip**: Use **absolute paths** for inputs.

`samples_tsv` formats:

- **FASTQ + paired-end** (`reads_type: "pe"`): `sample_id<TAB>R1<TAB>R2`
- **FASTQ + single-end** (`reads_type: "se"`): `sample_id<TAB>R1`
- **FASTQ + mixed** (`reads_type: "mx"`): `sample_id<TAB>R1<TAB>R2<TAB>U` (Bowtie2 mapping only)
- **BAM** (`format: "bam"`): `sample_id<TAB>bam_path`



## Outputs

Published outputs are written to `outdir` (default: `ChoCallate_output`), including standard Nextflow reports:

- `pipeline_report.html`
- `timeline_report.html`
- `trace.txt`

Consensus outputs depend on `output.type` and `output.format`:

- **Single merged** (default): `<outdir>/consensus.bcf` or `<outdir>/consensus.vcf.gz`
- **Per-sample**: `<outdir>/per_sample/<sample_id>.bcf` or `<sample_id>.vcf.gz`

Consensus VCF/BCF records include `FORMAT/GT`, `AD`, and `DP` (AD/DP from BAM), plus site-level `INFO` tags (`AN`, `AC`, `AF`, `NS`, total `DP`, etc.). Details: [Output Structure](https://github.com/alermol/ChoCallate/wiki/Output-Structure).

Consensus rule: every position of the coverage BED is written, monomorphic sites included. At each position the callers vote with their genotype (`REF`/`ALT`/`GT`), and the most frequent call wins when at least `consensus.threshold` callers support it. Positions without a caller record and hom-reference calls (any ploidy) vote for the reference allele; such positions are emitted as gVCF-style records with the symbolic ALT `<*>` (GT `0/0`, `INFO/END`, and `AD=(reference depth, other-allele depth)`), and caller records with a missing genotype (`./.`, `0/.`) abstain. Merged output uses `bcftools merge --gvcf`, which maps the `<*>` depth onto the real ALT at variant positions of other samples. `output.remove_invariant: true` drops the monomorphic records from the final output (per-sample files when `output.type: "sample"`, the merged file otherwise); with the default `false` all covered positions are kept.

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