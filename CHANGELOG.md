# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Changed

- Type of NM tag in resulting VCF was changed from A to R in order to count number of matching calls reference allele along with alternative ones.

### Added

- Add optional removing of monomorphic (invariant) positions both in single- and multi-sample output in order to reduce output file size.

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
