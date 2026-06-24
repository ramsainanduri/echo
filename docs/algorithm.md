# ECHO Algorithm

ECHO estimates CYP2D copy number and hybrid structure from per-base depth
profiles. It does not ingest BAM files and does not use CIGAR strings,
soft-clipping, split reads, or read-pair evidence.

## Inputs

The annotated target design BED is the coordinate reference. Its fourth column
is parsed dynamically to assign each interval to a gene and feature. For
example, labels such as `CYP2D6_NM_000106.6_ex9` and
`CYP2D7_NR_002570.6_ex9` define CYP2D6 and CYP2D7 exon targets without
hardcoded coordinates.

Depth input is separate. ECHO accepts samtools-style per-base depth
(`chrom pos depth`) and mosdepth-style per-base BED depth
(`chrom start end depth`).

## Tiling Normalization

Target panels can tile some genes more densely than others. ECHO takes a
gene-to-tiling file during PON construction:

```tsv
gene	tiling
CYP2D6	2
CYP2D7	2
```

For every target interval, ECHO divides the observed depth by the gene's tiling
factor. Genes not listed in the file are treated as 1X. The normalized target
baseline is then estimated from 1X adjusted targets.

## PON Construction

PON means panel of normals: a group of technically matched normal samples used
to estimate expected depth and variance across target intervals.

`build-pon` creates a panel-of-normals model:

1. Read normal samples from a manifest with `sample` and `depth_file` columns.
2. Summarize each per-base depth file over all annotated target intervals.
3. Apply gene tiling correction.
4. Normalize each sample to its adjusted 1X baseline.
5. Apply PCA correction to reduce systematic capture noise.
6. Store region-level PON mean and standard deviation.
7. Store PDS-level PON mean and standard deviation.
8. Infer sequencing modality from the variance of 1X background targets.

The PON stats output reports sample-level background depth, background CV,
covered target counts, PDS point counts, and serialized region/PDS statistics.

## BPSD Calling

`call-cnv` uses Bayesian Paralog Signal Deconvolution (BPSD). BPSD is ECHO's
depth-only approach for decomposing copy-number signal across homologous CYP2D
paralogs.

1. Normalize the test sample using the same target design and tiling map stored
   in the PON.
2. Estimate absolute copy number as `2 * sample_ratio / pon_ratio` for CYP
   regions.
3. Summarize calls at gene and exon levels.
4. Extract PDS depth signal. PDS means paralog-differentiating site: a position
   or target-base signal that helps separate homologous CYP2D paralogs. If the
   design does not contain explicit single-base PDS rows, ECHO falls back to
   dense CYP target-base depth signal.
5. Compare sample PDS depth to PON PDS depth to produce a PDS ratio and z-score.
6. Detect ratio shifts along genomic coordinate order.

## Bayesian Change-Point Detection

Change-point detection (CPD) identifies abrupt shifts in the ordered PDS ratio
signal. These shifts are reported as hybrid breakpoint candidates when they
have sufficient Bayesian support.

For an ordered PDS ratio signal `x`, ECHO compares:

- `M0`: one segment explains the whole signal
- `M1`: two segments split at candidate coordinate `t`

Within each segment, mean and variance are integrated out using a
Normal-Inverse-Gamma prior:

```text
sigma^2 ~ InvGamma(alpha0, beta0)
mu | sigma^2 ~ Normal(mu0, sigma^2 / kappa0)
```

The detector computes the log marginal likelihood for each segment and reports
the split with the strongest log Bayes factor (BF):

```text
log BF = log p(left) + log p(right) - log p(full)
```

Accepted splits are recursively refined, subject to minimum segment size and
log Bayes-factor thresholds.

In ECHO plots and call outputs, `BP` labels denote breakpoint candidates and
`BF` labels denote the Bayes-factor support for each candidate.
