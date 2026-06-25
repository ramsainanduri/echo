# ECHO

ECHO (Estimating Copy-number of Homologous Origins) is a pip-installable
Python package for per-base-depth CYP2D copy-number and hybrid breakpoint
calling using Bayesian Paralog Signal Deconvolution (BPSD).

ECHO uses an annotated target design BED to map target intervals to genes,
exons, introns, and CYP2D locus features. Depth profiles are provided
separately as per-base files from tools such as `samtools depth` or mosdepth.
The PON is built across every target in the BED. Calling is CYP2D-focused by
default, with optional copy-number summaries for additional standard genes
requested with `--genes`.

## Install

```bash
source .venv/bin/activate
pip install -e ".[dev]"
```

## Development

```bash
pre-commit install
ruff format .
ruff check .
mypy echo tests
pytest -q
```

Pull requests are expected to update `CHANGELOG.md` under `[Unreleased]`.

## Versioning

The package version is defined in [echo/version.py](echo/version.py):

```python
__version__ = "0.1.0"
```

Version bumps are manual and should happen in release pull requests. For normal
feature or bug-fix pull requests, add the user-facing change under
`[Unreleased]` in `CHANGELOG.md` and leave `echo/version.py` unchanged.

For a release pull request:

1. Update `echo/version.py`.
2. Move relevant changelog entries from `[Unreleased]` into a matching release
   heading, for example `## [0.2.0]`.
3. Verify the package metadata with `pip install -e .`.

> [!IMPORTANT]
> CI checks that release pull requests update `echo/version.py` and add a
> matching `CHANGELOG.md` release heading.

## CLI

```bash
echo-bpsd build-pon \
  --manifest examples/pon_samples.tsv \
  --bed examples/targets.bed \
  --tiling examples/tiling.tsv \
  --output examples/pon.pkl \
  --stats-output examples/pon.stats.json

echo-bpsd call-cnv \
  --depth examples/depths/sample.depth.bed \
  --pon examples/pon.pkl \
  --sample sample \
  --output-dir examples/ \
  --plot
```

Depth input may be samtools-style (`chrom pos depth`) or mosdepth-style
per-base BED (`chrom start end depth`).

> [!IMPORTANT]
> Depth files must use the same genome build and chromosome naming convention
> as the annotated BED.

The optional tiling file is a two-column TSV or CSV. Genes not listed are
treated as 1X:

```tsv
gene	tiling
CYP2D6	2
CYP2D7	2
```

> [!TIP]
> The tiling file should describe assay design, not expected biology. Genes
> absent from the file are assumed to be 1X tiled.

The PON sample manifest is a CSV or TSV with one row per normal:

```tsv
sample	depth_file
normal_001	depths/normal_001.depth.tsv
normal_002	depths/normal_002.depth.tsv
```

## Terminology

| Term | Meaning |
| --- | --- |
| CN | Copy number: the estimated number of genomic copies for a gene or exon. |
| CNV | Copy-number variation: deletion, duplication, or other dosage change. |
| PON | Panel of normals: technically matched normal samples used to model expected depth and variance. |
| PDS | Paralog-differentiating site: a position or target-base signal that helps distinguish homologous CYP2D paralogs. |
| DWBN | Density-Weighted Baseline Normalization: KDE-weighted baseline estimation for one-copy background targets. |
| BPSD | Bayesian Paralog Signal Deconvolution: ECHO's depth-only method for resolving paralog copy-number and hybrid structure. |
| CPD | Change-point detection: statistical detection of abrupt shifts in the ordered PDS ratio signal. |
| BF | Bayes factor: support for a two-segment breakpoint model over a one-segment model. |

### `build-pon` Arguments

| Argument | Required | Type / format | Description |
| --- | --- | --- | --- |
| `--manifest` | Required unless `--depth-dir` is used | CSV or TSV path | Table with `sample` and `depth_file` columns. This is the preferred input mode for reproducible PON builds and workflow engines. Relative depth paths are resolved relative to the manifest file. |
| `--depth-dir` | Required unless `--manifest` is used | Directory path | Directory containing normal depth files. Used as a convenience alternative to `--manifest`. |
| `--bed` | Yes | BED path | Annotated target design BED. The fourth column is parsed to assign intervals to genes and features. |
| `--tiling` | No | CSV or TSV path | Gene tiling-factor file with columns such as `gene` and `tiling`. Genes not listed are treated as 1X. |
| `--output` | Yes | `.pkl` path | Serialized PON model used by `call-cnv`. |
| `--stats-output` | No | `.json` or `.tsv` path | PON build QC and model summary. JSON contains full stats; TSV contains compact sample-level stats. |
| `--pattern` | No | Glob pattern | File pattern used with `--depth-dir`. Default: `*.tsv`. |
| `--max-pca-components` | No | Integer | Maximum number of PCA components used to reduce systematic capture noise. Default: `3`. |

### `call-cnv` Arguments

| Argument | Required | Type / format | Description |
| --- | --- | --- | --- |
| `--depth` | Yes | Depth TSV/BED path | Per-base depth profile for the sample being called. Supports `chrom pos depth` and `chrom start end depth`. |
| `--pon` | Yes | `.pkl` path | Serialized PON model produced by `build-pon`. Contains region statistics, PDS statistics, tiling factors, and modality metadata. |
| `--sample` | Yes | Text sample ID | Sample identifier used as the output filename prefix. |
| `--output-dir` | Yes | Directory path | Directory where ECHO writes the fixed output file set. Existing files with the same fixed names are overwritten. |
| `--genes` | No | Comma-separated gene symbols | Additional standard genes to call from the annotated BED, for example `TPMT,CYP2C19`. |
| `--plot` | No | Flag | Write CYP2D and standard-gene plots. When omitted, existing fixed plot files for the sample are removed. |
| `--plot-format` | No | `png`, `svg`, or `pdf` | Plot output format when `--plot` is used. Default: `png`. |

## Modes

### `build-pon`

Builds a serialized panel-of-normals model from normal samples. ECHO reads the
target design BED, summarizes each normal depth file over every design
interval, applies tiling correction, estimates each sample baseline with DWBN,
and stores PON statistics for all targets. CYP2D PDS statistics are stored
alongside the region statistics for hybrid breakpoint calling. The optional
`--stats-output` writes build QC metrics, including DWBN baseline values,
background coefficient of variation, covered regions, PDS points, and
serialized region/PDS statistics.

> [!IMPORTANT]
> Build the PON from technically comparable normal samples. Mixing capture
> designs, chemistry versions, or sequencing modalities can inflate variance
> and weaken CNV calls.

### `call-cnv`

Calls a single sample against a PON. ECHO normalizes sample depth using the
same target BED, tiling factors, and DWBN logic used during PON construction.
The primary call set is the CYP2D locus. Optional genes passed with `--genes`
are selected from the parsed BED and receive copy-number and z-score summaries
without CYP2D-specific PDS breakpoint analysis.

Outputs:

- `<sample>.echo.json`: full structured report
- `<sample>.cnv.tsv`: detailed gene, exon, and breakpoint calls
- `<sample>.genes.tsv`: compact gene-level copy-number table
- `<sample>.cyp2d6_cyp2d7.cnv.<format>`: CYP2D diagnostic plot when `--plot` is used
- `<sample>.<gene>.cnv.<format>`: one three-panel plot per requested standard gene when `--plot` is used

PON z-score is `(sample depth - PON mean) / PON SD`; negative values indicate
lower-than-expected depth and positive values indicate higher-than-expected
depth. Red vertical lines in the PDS panel are Bayesian change points. `BP`
means breakpoint candidate, and `BF` is the Bayes factor; larger BF values
indicate stronger support for a depth-ratio shift.

For standard genes present in the target BED and PON, pass a comma-separated
gene list:

```bash
echo-bpsd call-cnv \
  --depth sample.depth.bed \
  --pon pon.pkl \
  --sample sample \
  --output-dir results/ \
  --genes "TPMT,CYP2C19" \
  --plot \
  --plot-format svg
```

> [!TIP]
> Interpret CNV and breakpoint calls together with the diagnostic plot and QC
> metrics. Low coverage, sparse targets, or high PON variance can reduce
> confidence.

## Algorithm

ECHO operates on per-base depth, not BAM alignments. The high-level workflow is:

1. Parse the annotated BED to assign each target interval to a gene and feature.
2. Build or load PON statistics for every target interval in the BED.
3. Apply tiling correction and DWBN sample normalization.
4. Estimate target copy number as `2 * sample_normalized_depth / PON_mean`.
5. Summarize target calls to exon and gene copy number.
6. Extract PDS signal from explicit single-base CYP differentiating sites when
   present; otherwise use dense CYP target-base depth as a conservative depth
   signal.
7. Compare sample PDS ratios to PON PDS statistics.
8. Detect structural shifts with Bayesian change-point detection. Each
   candidate split compares a one-segment model to a two-segment model using a
   Normal-Inverse-Gamma marginal likelihood; accepted splits are reported as
   breakpoint candidates with log Bayes factors.

> [!NOTE]
> Explicit single-base PDS rows in the BED provide the most direct BPSD signal.
> If the design does not include them, ECHO falls back to dense CYP target-base
> depth, which is useful but less specific for paralog resolution.

> [!IMPORTANT]
> ECHO is research software. Clinical reporting requires local validation,
> quality thresholds, and review in the context of the full assay.

More details are in [docs/algorithm.md](docs/algorithm.md) and
[docs/outputs.md](docs/outputs.md).

## Runnable Example

Synthetic, non-sensitive inputs are provided in [examples/](examples/). The
example target BED is a compact subset of the real annotated design: a few
background regions plus the full CYP2D6/CYP2D7 block. The sample depth file
models a CYP2D6/2D7 hybrid-like shift with reduced downstream CYP2D7 depth.

```text
examples/
  targets.bed
  tiling.tsv
  pon_samples.tsv
  depths/
    normal_001.depth.bed
    normal_002.depth.bed
    sample.depth.bed
```

Run the example end to end:

```bash
echo-bpsd build-pon \
  --manifest examples/pon_samples.tsv \
  --bed examples/targets.bed \
  --tiling examples/tiling.tsv \
  --output examples/pon.pkl \
  --stats-output examples/pon.stats.json \
  --max-pca-components 0

echo-bpsd call-cnv \
  --depth examples/depths/sample.depth.bed \
  --pon examples/pon.pkl \
  --sample sample \
  --output-dir examples/ \
  --plot
```

The `.pkl`, `.json`, `.tsv`, and `.png` files are analysis outputs and can be
deleted or recreated at any time.

> [!NOTE]
> The runnable example uses only two normal samples, so PCA correction is
> disabled with `--max-pca-components 0`. Real PON builds should use a larger
> technically matched normal cohort before enabling PCA correction.
