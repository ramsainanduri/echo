# ECHO

ECHO (Estimating Copy-number of Homologous Origins) is a pip-installable
Python package for per-base-depth copy-number and CYP2D hybrid breakpoint
calling using Bayesian Paralog Signal Deconvolution (BPSD).

ECHO uses an annotated target design BED to map target intervals to genes,
exons, introns, and CYP2D locus features. Depth profiles are provided
separately as per-base files from tools such as `samtools depth` or mosdepth.

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
  --output examples/sample.echo.json \
  --cnv-output examples/sample.cnv.tsv \
  --gene-output examples/sample.genes.tsv \
  --plot examples/sample.cnv.png
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
| `--sample` | No | Text sample ID | Sample identifier used in the JSON report, plot title, and default output prefix. Defaults to the depth-file stem. |
| `--output` | No | `.json` or `.tsv` path | Main ECHO report. JSON is recommended for complete structured output; TSV is supported for simple tabular export. Defaults to `<sample>.echo.json`. |
| `--cnv-output` | No | TSV path | Dedicated CNV calls table with gene, exon, and breakpoint rows. |
| `--gene-output` | No | TSV path | Compact gene-level copy-number table with columns `gene`, `CN`, `CN(human)`, and `copy_number`. |
| `--plot` | No | `.png`, `.pdf`, or `.svg` path | Diagnostic CNV plot with absolute copy number, PON z-scores, PDS ratio signal, rolling median, and breakpoint annotations. |

## Modes

### `build-pon`

Builds a serialized panel-of-normals model from normal samples. ECHO reads the
target design BED, summarizes each normal depth file over the design intervals,
divides genes by their configured tiling factor, estimates the 1X baseline from
targets with 1X tiling, and stores PON region and PDS statistics. The optional
`--stats-output` writes build QC metrics, including sample-level background
depth, background coefficient of variation, covered regions, PDS points, and
serialized region/PDS statistics.

> [!IMPORTANT]
> Build the PON from technically comparable normal samples. Mixing capture
> designs, chemistry versions, or sequencing modalities can inflate variance
> and weaken CNV calls.

### `call-cnv`

Calls a single sample against a PON. ECHO normalizes sample depth using the
same tiling factors stored in the PON, estimates absolute copy number for CYP2D
regions, summarizes exon and gene copy number, and applies BPSD to detect PDS
ratio shifts consistent with hybrid breakpoints.

Outputs:

- `--output`: full JSON report, or TSV when the path ends in `.tsv`
- `--cnv-output`: dedicated TSV with gene, exon, and breakpoint calls
- `--gene-output`: compact gene-level copy-number TSV, for example
  `CYP2D6 3` and `CYP2D7 1`
- `--plot`: CNV diagnostic plot with a gene/exon model,
  absolute copy number, PON z-scores, PDS ratios, rolling median, and Bayesian
  breakpoint annotations

PON z-score is `(sample depth - PON mean) / PON SD`; negative values indicate
lower-than-expected depth and positive values indicate higher-than-expected
depth. Red vertical lines in the PDS panel are Bayesian change points. `BP`
means breakpoint candidate, and `BF` is the Bayes factor; larger BF values
indicate stronger support for a depth-ratio shift.

> [!TIP]
> Interpret CNV and breakpoint calls together with the diagnostic plot and QC
> metrics. Low coverage, sparse targets, or high PON variance can reduce
> confidence.

## Algorithm

BPSD operates on per-base depth, not BAM alignments. The workflow is:

1. Parse the annotated BED to assign each target interval to a gene and feature.
2. Extract PDS signal from explicit single-base CYP differentiating sites when
   present; otherwise use dense CYP target-base depth as a conservative depth
   signal.
3. Correct target depths by gene tiling factor, defaulting missing genes to 1X.
4. Normalize depths to the aggregate 1X baseline.
5. Compare sample PDS ratios to PON PDS statistics.
6. Detect structural shifts with Bayesian change-point detection. Each
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
  --output examples/sample.echo.json \
  --cnv-output examples/sample.cnv.tsv \
  --gene-output examples/sample.genes.tsv \
  --plot examples/sample.cnv.png
```

The `.pkl`, `.json`, `.tsv`, and `.png` files are analysis outputs and can be
deleted or recreated at any time.

> [!NOTE]
> The runnable example uses only two normal samples, so PCA correction is
> disabled with `--max-pca-components 0`. Real PON builds should use a larger
> technically matched normal cohort before enabling PCA correction.
