# Changelog

All notable changes to ECHO are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project uses semantic versioning.

## [Unreleased]

### Added

- Production Python package with `echo` and `echo-bpsd` command-line entry points.
- Panel-of-normals construction from sample manifests or depth-file directories.
- Density-Weighted Baseline Normalization (DWBN) for robust sample-level depth normalization.
- Optional tiling-factor correction for genes with non-uniform probe density.
- CYP2D6/CYP2D7-focused copy-number calling with exon-level, gene-level, and hybrid breakpoint outputs.
- Bayesian Paralog Signal Deconvolution for CYP2D6/CYP2D7 hybrid detection.
- Optional standard-gene copy-number calling through `call-cnv --genes`.
- Standardized output layout with structured JSON reports, detailed CNV TSV files, compact gene CN TSV files, PON statistics, and diagnostic plots.
- Publication-style CYP2D and standard-gene diagnostic plots with PNG, SVG, and PDF output support.
- Continuous integration checks for formatting, linting, typing, tests, changelog updates, and release version consistency.

### Documentation

- Added usage documentation, argument tables, output schemas, and algorithm details for PON construction, DWBN normalization, CYP2D CNV calling, optional standard-gene calling, and plotting.
