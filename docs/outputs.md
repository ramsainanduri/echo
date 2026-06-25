# ECHO Outputs

## PON Outputs

`build-pon --output pon.pkl` writes the serialized model used by `call-cnv`.
The PON is built for every target in the annotated BED, not only CYP2D targets.
The model contains:

- parsed target design regions
- tiling factors
- region-level PON means and standard deviations for all BED targets
- CYP2D PDS-level PON means and standard deviations
- PCA correction components
- inferred sequencing modality
- sample-level build statistics

`build-pon --stats-output pon.stats.json` writes build statistics. When the
path ends in `.tsv`, ECHO writes the compact sample-level stats table.

Sample-level DWBN fields include:

- `dwbn_1x_baseline_depth`
- `dwbn_kde_mode_depth`
- `dwbn_effective_region_count`
- `dwbn_min_weight`
- `dwbn_max_weight`
- `background_raw_regions`
- `background_weighted_regions`
- `background_cv`

## Call Outputs

`call-cnv --sample sample --output-dir results/` writes a fixed text output
set using `sample` as the filename prefix. Plot files are written only when
`--plot` is used.

## Call Report

`results/sample.echo.json` contains the complete structured report:

- sample and PON metadata
- tiling factors
- optional `standard_genes` requested with `--genes`
- gene-level copy-number calls
- exon or target-level copy-number calls
- CYP2D breakpoint candidates
- CYP2D hybrid candidate rows
- quality metrics

## CNV Calls TSV

`results/sample.cnv.tsv` contains a flat table with stable columns:

```text
call_type
gene
copy_number
integer_copy_number
hybrid
feature
chrom
start
end
z_score
breakpoint_coordinate
log_bayes_factor
left_pds_ratio
right_pds_ratio
segments
```

Rows use `call_type=gene`, `call_type=exon`, or `call_type=hybrid`.
Hybrid rows are CYP2D-specific and include breakpoint support fields.

## Gene Summary TSV

`results/sample.genes.tsv` contains a compact gene-level table:

```tsv
gene	CN	copy_number
CYP2D6	3	2.94
CYP2D7	1	1.12
```

When `--genes` is used, requested standard genes are included in the same file
if they are present in the PON BED and have callable target intervals.

## Plots

`results/sample.cyp2d6_cyp2d7.cnv.<format>` contains the CYP2D diagnostic plot
when `--plot` is used. The plot format is selected with `--plot-format` and
defaults to `png`:

1. CYP2D gene/exon model.
2. Absolute copy-number estimates by target interval.
3. PON z-score scatter with z-score reference lines.
4. PDS ratio signal and Bayesian breakpoint annotations.

`call-cnv --genes TPMT,CYP2C19 --plot` writes one standard gene plot per
requested gene in the same output directory, for example
`results/sample.tpmt.cnv.png` or `results/sample.tpmt.cnv.svg`. Each
standard-gene plot has three panels:

1. Gene model.
2. Absolute copy number by target interval.
3. PON z-score by target interval.

PDS and hybrid breakpoint panels are only generated for the CYP2D locus.
When `--plot` is omitted, ECHO removes existing fixed plot files for that
sample so stale plot files are not left beside fresh text outputs.
