# ECHO Algorithm

ECHO is a CYP2D-focused copy-number caller for targeted panels and similar
per-base depth assays. The panel of normals (PON) is built across every region
in the annotated target BED so that normalization and variance modeling use the
full assay design. During sample calling, ECHO always evaluates the CYP2D locus
as the primary target. Additional genes can be requested with `--genes`; those
standard genes are selected dynamically from the same BED and are reported with
copy-number and z-score summaries, but they do not use CYP2D-specific PDS
breakpoint logic.

## Inputs

The annotated target design BED is the coordinate and annotation source. ECHO
parses the fourth BED column to derive gene and feature labels, so labels such
as `CYP2D6_NM_000106.6_ex9`, `CYP2D7_NR_002570.6_ex9`, or
`TPMT_NM_000367.5_ex4` determine which targets belong to each gene. The PON and
caller use the parsed BED regions; gene coordinates are not supplied in code.

Depth input is separate from the BED. ECHO accepts both samtools-style depth
(`chrom pos depth`) and BED-like per-base depth (`chrom start end depth`). The
BED and depth files must use the same genome build and chromosome naming.

## Tiling Correction

Target panels can tile genes at different densities. During PON construction,
ECHO can read a gene tiling table:

```tsv
gene	tiling
CYP2D6	2
CYP2D7	2
```

For every target interval, observed depth is divided by that gene's tiling
factor. Genes absent from the tiling table are treated as 1X. This correction is
applied before sample-level baseline normalization.

## DWBN Baseline Normalization

ECHO uses Density-Weighted Baseline Normalization (DWBN) to estimate each
sample's one-copy background depth. DWBN reduces the influence of off-center
background targets without hard exclusion thresholds.

For each sample:

1. Summarize mean depth over every target interval in the BED.
2. Apply gene tiling correction to obtain adjusted target depths.
3. Select 1X background targets. If no explicit 1X background targets are
   present, use all adjusted targets as the fallback baseline pool.
4. Fit a Gaussian KDE to the adjusted background depths.
5. Evaluate the KDE over the observed background depths and use those density
   values as trust weights.
6. Calculate the sample baseline as the density-weighted average of adjusted
   background depths.
7. Normalize every target interval by this DWBN baseline.

The KDE mode is recorded as a QC value, but the baseline used for normalization
is the KDE-density-weighted average. The effective weighted region count helps
show whether the baseline was supported by many similarly trusted targets or by
a narrow subset of the background pool.

## PON Construction

`build-pon` creates a serialized model for the full target BED:

1. Read normal samples from a manifest with `sample` and `depth_file` columns,
   or from a depth directory.
2. Parse all BED regions and assign each interval to a gene and feature.
3. Summarize each normal depth file over all BED regions.
4. Apply tiling correction and DWBN normalization per sample.
5. Optionally apply PCA correction to reduce systematic capture noise.
6. Store region-level PON means and standard deviations for every target.
7. Extract and normalize CYP2D PDS signal and store PDS-level PON statistics.
8. Record sample-level QC fields, including DWBN baseline depth, KDE mode,
   effective weighted background count, background CV, covered target count,
   and PDS point count.

Because the PON stores statistics for all target regions, the same PON can be
used later to call CYP2D and any optional standard genes present in the BED.

## CNV Calling

`call-cnv` normalizes one test sample using the target BED, tiling map, and PON
statistics stored in the model.

For each called target interval:

```text
sample_ratio = sample_normalized_depth / PON_mean_normalized_depth
copy_number = 2 * sample_ratio
z_score = (sample_normalized_depth - PON_mean_normalized_depth) / PON_sd
```

Copy number is summarized at two levels:

- target/exon level: per-interval copy number and PON z-score
- gene level: median target copy number, rounded integer CN, and target count

CYP2D6, CYP2D7, and CYP2D8P are always selected as the core locus. Genes passed
with `--genes` are added to the selected target set if their parsed gene symbols
exist in the PON BED. Missing requested genes fail fast with a clear error.

## CYP2D Hybrid Detection

CYP2D calling includes Bayesian Paralog Signal Deconvolution (BPSD). BPSD uses
PDS signal to detect shifts consistent with CYP2D6/CYP2D7 hybrid structure.
PDS means paralog-differentiating site: a position or target-base signal that
helps separate homologous CYP2D paralogs. If explicit single-base PDS rows are
not present in the BED, ECHO falls back to dense CYP2D target-base depth.

For the ordered PDS ratio signal, ECHO compares:

- `M0`: one segment explains the full signal
- `M1`: two segments split at a candidate coordinate

Within each segment, mean and variance are integrated out using a
Normal-Inverse-Gamma prior:

```text
sigma^2 ~ InvGamma(alpha0, beta0)
mu | sigma^2 ~ Normal(mu0, sigma^2 / kappa0)
```

The detector reports candidate split points by log Bayes factor:

```text
log BF = log p(left) + log p(right) - log p(full)
```

Accepted splits are recursively refined subject to minimum segment size and log
Bayes-factor thresholds. Standard genes requested with `--genes` do not use
PDS or BPSD; they receive copy-number and PON z-score calls only.

## Plots

The CYP2D plot has four panels: gene model, absolute copy number, PON z-score,
and PDS/BPSD signal. Optional standard genes receive one three-panel PNG per
gene: gene model, absolute copy number, and PON z-score.
