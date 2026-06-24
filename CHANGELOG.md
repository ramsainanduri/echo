# Changelog

All notable changes to ECHO are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project uses semantic versioning.

## [Unreleased]

### Added

- Initial pip-installable ECHO package and `echo-bpsd` command-line interface.
- Panel-of-normals construction from sample manifests and per-base depth files.
- Tiling-aware depth normalization using configurable gene tiling factors.
- CYP2D copy-number calling with exon-level reports and breakpoint detection.
- Bayesian Paralog Signal Deconvolution with Bayesian change-point detection.
- Optional CNV diagnostic plots with copy-number, z-score, and PDS-ratio panels.
- Dedicated CNV calls TSV output and PON build statistics output.
- Compact gene-level copy-number TSV output.
- CNV output with hybrid candidate rows.
- CNV diagnostic plots with gene/exon model, ordered exon axis, z-score color
  scale, and breakpoint annotations.
- `call-cnv --sample` for stable sample IDs, plot titles, and default output
  filename prefixes.
- Expanded acronym documentation.
- Algorithm and output documentation for PON construction and CNV calling.
- Development quality checks for formatting, linting, typing, tests, and
  changelog enforcement.
- Release version consistency checks in GitHub Actions.
