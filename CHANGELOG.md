# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

## [3.2.3] - 2026-10-01

### Fixed

- Each sample's BAM is now joined with its own coverage BED by `sample_id` before variant calling and consensus generation. `GENERATE_COVERAGE` no longer drops `sample_id` from its output, so the coverage BED can no longer be paired with another sample's BAM and callers no longer report variants at positions below `coverage.min_coverage`.
- Merged multi-sample outputs now have a deterministic sample order that follows `samples_tsv`. `MERGE_OUTPUTS` ranks the per-sample consensus files by their `samples_tsv` position instead of relying on process completion order, and the merge list is sorted naturally (`ls -1v`) so the order is preserved for ≥ 10 samples.
- BAM input (`format: "bam"`) now works: the sample channel no longer reads a separate `input.input_format` key that disagreed with the `input.format` key used everywhere else.
- `CLIParamsValidation.effective_callers_validation` now actually validates caller names and no longer calls `split()` before its null check; `cons_threshold_validation` also checks `callers` for null before using it.

### Changed

- The documented input format key is now `format` (was `input_format`) in the README, matching `assets/templates/config.yaml`.

## [3.2.2] - 2026-09-30

### Fixed

- `gatk SplitIntervals` in all calling processes (`CALLING_BCFTOOLS`, `CALLING_FREEBAYES`, `CALLING_GATK`) now writes its temp files to the process `$TMPDIR` (`--java-options "-Djava.io.tmpdir=$TMPDIR"` and `--tmp-dir $TMPDIR`) instead of the system temporary directory.

## [3.2.1] - 2026-09-30

### Changed

- The `bcftools` caller now emits `FORMAT/AD` and `FORMAT/DP`, which are requested directly from `bcftools mpileup` (`-a FORMAT/AD,FORMAT/DP`) instead of being recovered downstream.
- The reference sequence dictionary (`ref_genome.dict`) is now passed to all calling processes and to `GENERATE_CONSENSUS`, as required by `gatk SplitIntervals`.
- All calling processes (`CALLING_BCFTOOLS`, `CALLING_FREEBAYES`, `CALLING_GATK`) and `GENERATE_CONSENSUS` now split the coverage BED with `gatk SplitIntervals` (balanced, size-aware intervals) instead of `split`/`bedops --chop` + `split`.
- `FORMAT/AD` and `FORMAT/DP` recovery in `generate_consensus.py` is faster: the genotype walk resumes from the last record instead of rescanning from the start at every position, and the BAM pileup is streamed position-by-position instead of materialising the whole region (bounded memory, no second pass).

### Fixed

- Removal of monomorphic (invariant) positions now filters hom-ref sites with `COUNT(GT="RR")=N_SAMPLES` instead of counting all-homozygous/all-heterozygous genotypes, so only truly invariant positions are dropped and genuine variants are no longer removed.

## [3.2.0] - 2026-09-04

### Added

- Consensus outputs now include `FORMAT/AD` and `FORMAT/DP`, recovered from the sample BAM pileup (overlapping PE mates are deduplicated for allele counts; `DP` counts all covering reads to match `samtools depth` used for the coverage BED).
- Site-level `INFO` annotations via `bcftools +fill-tags` on both per-sample and merged outputs: `AN`, `AC`, `AF`, `NS`, `AC_Hom`, `AC_Het`, `MAF`, `TYPE`, `F_MISSING`, and `DP` (sum of `FORMAT/DP`).

### Changed

- Consensus `QUAL` is always missing (`.`) instead of a placeholder score.
- `GENERATE_CONSENSUS` takes the sample BAM as input so AD/DP can be recovered during consensus writing.

## [3.1.0] - 2026-09-01

### Added

- Wiki: [Merge Strategies and Testing](https://github.com/alermol/ChoCallate/wiki/Merge-Strategies-and-Testing).

### Changed

- `MERGE_OUTPUTS` now uses chunked pipe merge for N samples ≥ 2 (no intermediate files, named pipes only).

## [3.0.0] - 2026-08-27

### Changed

- Dependency version pins in `environment.yaml` were relaxed to minimum versions (`>=`), and the Nextflow requirement was raised to `>=26.04.6`.
- The merging of single-sample VCF/BCF files into multi-sample files was significantly accelerated by using hierarchical merging.
- `MERGE_OUTPUTS` now uses `consensus.cpu` for thread allocation. The previous retry strategy that decreased CPU count on each attempt was removed.
- Default `output.type` changed from `"sample"` to `"single"` (one multi-sample consensus file).

### Added

- Optional removal of monomorphic (invariant) positions in both single- and multi-sample outputs to reduce file size (`output.remove_invariant`).
- Optional splitting of multi-allelic variants in output multi-sample VCF/BCF files (`output.split_multiallelic`).
- Added information to the output VCF/BCF file headers about the ChoCallate version (`##tool`), multi-allelic variant splitting (`##multialleleSplit`), and monomorphic variant removal (`##noMonomorphic`). The consensus threshold header key was renamed to `##consensusThreshold`.

### Removed

- The NM tag was temporarily removed due to issues during merging. It will be restored in future updates. The consensus threshold remains in the VCF/BCF header.

### Fixed

- Fixed a bug that caused GATK to crash when the system's temporary directory became full.

## [2.0.5] - 2026-05-20

### Fixed

- Fixed typos in the MAPPING_BWA and MAPPING_MINIMAP2 processes that caused them to crash when launched.

### Added

- The MERGE_OUTPUTS process now uses a retry strategy with a decreasing number of processors for each attempt.

## [2.0.4] - 2026-05-09

### Changed

- Refractored the pipeline code and file structure towards nf-core standards.

### Fixed

- Fixed hardcoded paths to input files in the config for test run.

## [2.0.3] - 2026-04-27

### Changed

- Removed `maxForks 1` in several steps, including BAM file filtering, BED file generation with coverage depth, duplicate removal, left-alignment of indels, and BCF to VCF conversion. These were performance bottlenecks for running on many samples.
- Cleanup observer for nf-boost has been changed to `'v2'`.

## [2.0.2] - 2026-04-25

### Added

- Implemented launch in a Docker container.

## [2.0.1] - 2026-04-22

### Fixed

- Name collision was fixed in the coverage generation process, arising from using the same `NO_FILE` for both include_bed and exclude_bed.
- During large file processing, the system temporary directory could overflow, which would result in the pipeline crashing. For each process, the system `$TMPDIR` variable is overridden with a temporary directory that uses the process's `$PWD`.

## [2.0.0] - 2026-04-21

### Changed

- Sequential BCF/VCF indexing was replaced with parallel indexing when merging output files for each sample into a single multi-sample file
- The version number was changed to 2.0.0 due to the removal of backward compatibility in version 1.1.0 compared to 1.0.1, in accordance with semantic versioning

## [1.1.0] - 2026-04-19

### Added

- The `include bed` and `exclude bed` arguments have been added as a more flexible replacement for the `custom bed` argument

### Changed

- All test run data has been moved to a separate `test_run` directory to clean up the project root folder

### Removed

- Argument `custom_bed` has been removed

## [1.0.1] - 2026-04-16

### Fixed

- Added the missing validation of the `reference_index_dir` parameter, which was declared as required, but its absence caused a failure during reads mapping instead of an early termination of the pipeline

## [1.0.0] - 2026-04-14

Initial release
