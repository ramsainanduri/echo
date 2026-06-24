# ECHO Outputs

## PON Outputs

`build-pon --output pon.pkl` writes the serialized model used by `call-cnv`.
The model contains:

- parsed target design regions
- tiling factors
- region-level PON means and standard deviations
- PDS-level PON means and standard deviations
- PCA components
- inferred sequencing modality
- sample-level build statistics

`build-pon --stats-output pon.stats.json` writes a human-readable build stats
file. When the output path ends in `.tsv`, ECHO writes the compact sample-level
stats table.

## Call Outputs

`call-cnv --output sample.echo.json` writes the full report:

- sample and PON metadata
- tiling factors
- gene copy-number calls
- exon copy-number calls
- breakpoint candidates
- quality metrics

`call-cnv --cnv-output sample.cnv.tsv` writes a flat calls table with `level`
values of `gene`, `exon`, or `breakpoint`.

`call-cnv --plot sample.cnv.png` writes a diagnostic figure with three panels:

1. Absolute copy-number estimates over genomic coordinates.
2. PON z-score scatter with threshold guide lines.
3. PDS ratio signal, rolling median, and Bayesian breakpoint annotations.

The plot is designed for scientific review of the call: it shows the copy
number estimate, whether the call is supported by PON-normalized deviation, and
where the PDS signal changes along the CYP locus.
