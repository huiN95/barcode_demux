# Changelog

All notable changes to this project will be documented in this file.

## [1.0.5] - 2026-04-28

### Changed
- Updated version to 1.0.5 across the project.
- Updated Conda recipe for version 1.0.5.

## [1.0.4] - 2026-04-28

### Added
- **Flexible CLI Arguments**: Added underscore aliases for all long command-line arguments (e.g., `--input_file` as an alias for `--input-file`).
- Updated documentation for version 1.0.4.

### Changed
- Updated version to 1.0.4 in `Cargo.toml` and `src/cli.rs`.

## [1.0.3] - 2026-04-28

### Added
- Pipeline version 4 support with `min-pair-len` and `min-pair-score`.
- `keep-barcode` and `single-end-filter` options.

## [1.0.2] - 2026-04-28

### Added
- Comprehensive `README.md` documentation.
- Standard MIT License and Research-only disclaimer.

## [1.0.1] - Previous Version

### Added
- Initial support for BAM format.
- Myers bit-parallel algorithm integration.

## [1.0.0] - Initial Release

- Core demultiplexing pipeline.
- Support for FASTA/FASTQ formats.
