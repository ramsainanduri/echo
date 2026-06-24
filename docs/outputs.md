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

`build-pon --stats-output pon.stats.json` writes build statistics. When the
output path ends in `.tsv`, ECHO writes the compact sample-level stats table.

## Call Outputs

`call-cnv --output sample.echo.json` writes the full report:

- sample and PON metadata
- tiling factors
- gene copy-number calls
- exon copy-number calls
- breakpoint candidates
- quality metrics

`call-cnv --cnv-output sample.cnv.tsv` writes a flat calls table. The first
columns are `call_type`, `gene`, and `CN(human)`. Hybrid breakpoint candidates
are emitted as `call_type=hybrid` rows with breakpoint coordinates and
Bayes-factor support.

`call-cnv --gene-output sample.genes.tsv` writes the compact gene-level copy
number table:

```tsv
gene	CN	CN(human)	copy_number
CYP2D6	3	3 copies (estimated 2.94)	2.94
CYP2D7	1	1 copy (estimated 1.12)	1.12
```

`call-cnv --plot sample.cnv.png` writes a diagnostic figure with four panels:

1. CYP2D gene/exon model.
2. Absolute copy-number estimates by target interval.
3. PON z-score scatter with threshold guide lines and a z-score color scale.
4. PDS ratio signal, rolling median, ratio guide lines, and the strongest
   Bayesian breakpoint annotations.

The plot is designed for scientific review of the call: it shows the copy
number estimate, whether the call is supported by PON-normalized deviation, and
where the PDS signal changes along the CYP locus. The x-axis is labeled by CYP
gene/exon positions rather than raw coordinates where exon annotations are
available.

PON z-score is `(sample depth - PON mean) / PON SD`. Values near zero indicate
that the sample resembles the panel of normals at that target; strongly
negative values support lower copy number, and strongly positive values support
higher copy number.

In the PDS panel, red vertical lines are Bayesian change points. Labels such as
`BP1 BF 7179` mean breakpoint candidate 1 with Bayes-factor support of 7179.
Higher BF values indicate stronger support for a two-segment depth-ratio model
than for a single unchanged segment.
