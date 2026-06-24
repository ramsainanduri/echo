"""Depth-file loading and interval summarisation."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd

from echo.regions import BedRegion


class DepthProfile:
    """Per-base depth profile indexed by chromosome and zero-based position."""

    def __init__(self, frame: pd.DataFrame, source: str | Path) -> None:
        self.frame = frame
        self.source = str(source)

    @classmethod
    def from_path(cls, path: str | Path) -> DepthProfile:
        """Load a samtools-depth or mosdepth per-base depth file.

        Parameters
        ----------
        path
            Input TSV. Supported shapes are ``chrom pos depth`` with one-based
            positions and ``chrom start end depth`` with zero-based BED
            intervals.

        Returns
        -------
        DepthProfile
            Normalized per-base depth table with columns ``chrom``, ``pos0``,
            and ``depth``.
        """

        input_path = Path(path)
        raw = pd.read_csv(input_path, sep="\t", header=None, comment="#")
        if raw.shape[1] < 3:
            raise ValueError(f"Depth file {path} must contain at least three columns")
        if raw.shape[1] >= 4:
            frame = cls._from_bed_depth(raw)
        else:
            frame = pd.DataFrame(
                {
                    "chrom": raw.iloc[:, 0].astype(str),
                    "pos0": raw.iloc[:, 1].astype(np.int64) - 1,
                    "depth": raw.iloc[:, 2].astype(float),
                }
            )
        frame = frame.sort_values(["chrom", "pos0"]).reset_index(drop=True)
        return cls(frame, input_path)

    @staticmethod
    def _from_bed_depth(raw: pd.DataFrame) -> pd.DataFrame:
        """Expand or compact mosdepth BED rows to per-base positions."""

        pieces: list[pd.DataFrame] = []
        for chrom, start, end, depth in raw.iloc[:, :4].itertuples(index=False, name=None):
            start_i = int(start)
            end_i = int(end)
            if end_i <= start_i:
                continue
            positions = np.arange(start_i, end_i, dtype=np.int64)
            pieces.append(
                pd.DataFrame(
                    {
                        "chrom": str(chrom),
                        "pos0": positions,
                        "depth": float(depth),
                    }
                )
            )
        if not pieces:
            raise ValueError("Depth BED contained no positive-length intervals")
        return pd.concat(pieces, ignore_index=True)

    def depths_for_region(self, region: BedRegion) -> np.ndarray:
        """Return all observed per-base depths within a BED region."""

        subset = self.frame[
            (self.frame["chrom"] == region.chrom)
            & (self.frame["pos0"] >= region.start)
            & (self.frame["pos0"] < region.end)
        ]
        return subset["depth"].to_numpy(dtype=float)

    def mean_depth(self, region: BedRegion) -> float:
        """Return mean depth over one interval or NaN when uncovered."""

        depths = self.depths_for_region(region)
        if depths.size == 0:
            return float("nan")
        return float(np.nanmean(depths))

    def region_means(self, regions: Iterable[BedRegion]) -> pd.DataFrame:
        """Compute mean depth for each interval."""

        rows = []
        for region in regions:
            rows.append(
                {
                    "chrom": region.chrom,
                    "start": region.start,
                    "end": region.end,
                    "name": region.name,
                    "gene": region.gene,
                    "feature": region.feature,
                    "mean_depth": self.mean_depth(region),
                }
            )
        return pd.DataFrame(rows)

    def pds_signal(
        self, regions: list[BedRegion], fallback_regions: list[BedRegion]
    ) -> pd.DataFrame:
        """Extract PDS depth signal, falling back to dense CYP target bases.

        Parameters
        ----------
        regions
            Explicit single-base PDS regions.
        fallback_regions
            CYP target regions used when the design has no explicit PDS rows.

        Returns
        -------
        pandas.DataFrame
            Columns ``chrom``, ``pos0``, ``depth``, and ``source``.
        """

        if regions:
            rows = []
            for region in regions:
                depths = self.depths_for_region(region)
                if depths.size:
                    rows.append(
                        {
                            "chrom": region.chrom,
                            "pos0": region.start,
                            "depth": float(np.nanmean(depths)),
                            "source": region.name,
                        }
                    )
            return pd.DataFrame(rows).sort_values(["chrom", "pos0"]).reset_index(drop=True)

        mask = pd.Series(False, index=self.frame.index)
        for region in fallback_regions:
            mask |= (
                (self.frame["chrom"] == region.chrom)
                & (self.frame["pos0"] >= region.start)
                & (self.frame["pos0"] < region.end)
            )
        result = self.frame.loc[mask, ["chrom", "pos0", "depth"]].copy()
        result["source"] = "dense_cyp_target_fallback"
        return result.sort_values(["chrom", "pos0"]).reset_index(drop=True)
