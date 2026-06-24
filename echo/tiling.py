"""Gene tiling-factor parsing and normalization helpers."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd

from echo.regions import BedRegion


def load_tiling_factors(path: str | Path | None) -> dict[str, float]:
    """Load gene tiling factors from a two-column text file.

    The file may be tab- or comma-delimited, with or without a header. Accepted
    header names are ``gene`` and one of ``tiling``, ``factor``, or
    ``tiling_factor``. Genes not present in the returned mapping are treated as
    1X by the normalization helpers.

    Parameters
    ----------
    path
        Optional tiling file path.

    Returns
    -------
    dict
        Mapping from gene symbol to positive tiling factor.
    """

    if path is None:
        return {}
    tiling_path = Path(path)
    with tiling_path.open("r", encoding="utf-8", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
        rows = list(csv.reader(handle, dialect))
    rows = [row for row in rows if row and not row[0].startswith("#")]
    if not rows:
        return {}

    first = [cell.strip().lower() for cell in rows[0]]
    has_header = "gene" in first and any(
        name in first for name in ("tiling", "factor", "tiling_factor")
    )
    data_rows = rows[1:] if has_header else rows
    if has_header:
        gene_idx = first.index("gene")
        factor_idx = next(
            first.index(name) for name in ("tiling", "factor", "tiling_factor") if name in first
        )
    else:
        gene_idx = 0
        factor_idx = 1

    factors: dict[str, float] = {}
    for line_number, row in enumerate(data_rows, start=2 if has_header else 1):
        if len(row) <= max(gene_idx, factor_idx):
            raise ValueError(f"Tiling file line {line_number} has fewer than two columns")
        gene = row[gene_idx].strip()
        factor = float(row[factor_idx])
        if not gene:
            raise ValueError(f"Tiling file line {line_number} has an empty gene name")
        if not np.isfinite(factor) or factor <= 0.0:
            raise ValueError(f"Tiling factor for {gene!r} on line {line_number} must be positive")
        factors[gene] = factor
    return factors


def tiling_factor_for_gene(gene: str, tiling_factors: dict[str, float]) -> float:
    """Return a gene tiling factor, defaulting to one."""

    return float(tiling_factors.get(gene, 1.0))


def add_region_tiling(region_means: pd.DataFrame, tiling_factors: dict[str, float]) -> pd.DataFrame:
    """Attach tiling factors and adjusted depths to a region-depth dataframe."""

    adjusted = region_means.copy()
    adjusted["tiling_factor"] = adjusted["gene"].map(
        lambda gene: tiling_factor_for_gene(str(gene), tiling_factors)
    )
    adjusted["adjusted_depth"] = adjusted["mean_depth"] / adjusted["tiling_factor"]
    return adjusted


def baseline_from_adjusted_regions(adjusted_regions: pd.DataFrame) -> float:
    """Estimate the one-copy baseline from adjusted target depths."""

    background = adjusted_regions.loc[
        adjusted_regions["tiling_factor"] == 1.0, "adjusted_depth"
    ].dropna()
    if background.empty:
        background = adjusted_regions["adjusted_depth"].dropna()
    if background.empty:
        raise ValueError("No covered target regions are available for baseline normalization")
    return float(np.nanmedian(background))


def normalize_region_means(
    region_means: pd.DataFrame, tiling_factors: dict[str, float]
) -> pd.Series:
    """Normalize target depths after correcting each gene for tiling factor."""

    adjusted = add_region_tiling(region_means, tiling_factors)
    baseline = baseline_from_adjusted_regions(adjusted)
    values = adjusted["adjusted_depth"] / max(baseline, np.finfo(float).eps)
    return pd.Series(values.to_numpy(dtype=float), index=adjusted["name"].tolist())


def pds_tiling_factor(
    chrom: str, pos0: int, regions: list[BedRegion], tiling_factors: dict[str, float]
) -> float:
    """Return the tiling factor for a PDS coordinate."""

    for region in regions:
        if region.chrom == chrom and region.start <= pos0 < region.end:
            return tiling_factor_for_gene(region.gene, tiling_factors)
    return 1.0


def normalize_pds_signal(
    pds: pd.DataFrame,
    region_means: pd.DataFrame,
    regions: list[BedRegion],
    tiling_factors: dict[str, float],
) -> pd.DataFrame:
    """Normalize PDS depths by gene tiling and one-copy target baseline."""

    if pds.empty:
        return pds
    adjusted_regions = add_region_tiling(region_means, tiling_factors)
    baseline = baseline_from_adjusted_regions(adjusted_regions)
    normalized = pds.copy()
    normalized["tiling_factor"] = [
        pds_tiling_factor(str(row.chrom), int(cast(Any, row.pos0)), regions, tiling_factors)
        for row in normalized.itertuples(index=False)
    ]
    normalized["depth"] = normalized["depth"] / normalized["tiling_factor"]
    normalized["depth"] = normalized["depth"] / max(baseline, np.finfo(float).eps)
    return normalized


def serializable_tiling_factors(tiling_factors: dict[str, float]) -> dict[str, Any]:
    """Return a JSON/pickle-friendly sorted tiling map."""

    return {gene: float(tiling_factors[gene]) for gene in sorted(tiling_factors)}
