"""Panel-of-normals model construction and serialization."""

from __future__ import annotations

import json
import pickle
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from echo.depth import DepthProfile
from echo.manifest import DepthSample
from echo.regions import AnnotatedBed, BedRegion
from echo.tiling import (
    add_region_tiling,
    load_tiling_factors,
    normalize_pds_signal,
    normalize_region_means,
    serializable_tiling_factors,
)
from echo.version import __version__


@dataclass
class PONModel:
    """Serialized ECHO panel-of-normals model."""

    version: str
    bed_regions: list[BedRegion]
    region_stats: pd.DataFrame
    pds_stats: pd.DataFrame
    modality: str
    pca_components: np.ndarray
    region_order: list[str]
    pds_mode: str
    tiling_factors: dict[str, float]
    sample_stats: pd.DataFrame

    def save(self, path: str | Path) -> None:
        """Serialize the model with pickle."""

        with Path(path).open("wb") as handle:
            pickle.dump(self, handle, protocol=pickle.HIGHEST_PROTOCOL)

    @classmethod
    def load(cls, path: str | Path) -> PONModel:
        """Load a serialized model."""

        with Path(path).open("rb") as handle:
            loaded = pickle.load(handle)
        if not isinstance(loaded, cls):
            raise TypeError(f"{path} does not contain a PONModel")
        return loaded

    @property
    def bed(self) -> AnnotatedBed:
        """Return the parsed BED wrapper."""

        return AnnotatedBed(self.bed_regions)

    def build_stats(self) -> dict[str, object]:
        """Return JSON-serializable PON build statistics."""

        return {
            "version": self.version,
            "modality": self.modality,
            "pds_mode": self.pds_mode,
            "sample_count": int(len(self.sample_stats)),
            "target_region_count": int(len(self.region_stats)),
            "pds_site_count": int(len(self.pds_stats)),
            "pca_component_count": int(self.pca_components.shape[0]),
            "tiling_factors": self.tiling_factors,
            "sample_stats": self.sample_stats.to_dict(orient="records"),
            "region_stats": self.region_stats.to_dict(orient="records"),
            "pds_stats": self.pds_stats.to_dict(orient="records"),
        }

    def write_stats(self, path: str | Path) -> None:
        """Write PON build statistics as JSON or a compact TSV summary."""

        out = Path(path)
        if out.suffix.lower() == ".tsv":
            self.sample_stats.to_csv(out, sep="\t", index=False)
            return
        with out.open("w", encoding="utf-8") as handle:
            json.dump(self.build_stats(), handle, indent=2, sort_keys=True)
            handle.write("\n")


class PanelOfNormalsBuilder:
    """Build a PCA-corrected panel of normals from per-base depths."""

    def __init__(
        self,
        bed_path: str | Path,
        tiling_path: str | Path | None = None,
        max_pca_components: int = 3,
        targeted_variance_cv: float = 0.18,
    ) -> None:
        self.bed = AnnotatedBed.from_path(bed_path)
        self.tiling_factors = load_tiling_factors(tiling_path)
        self.max_pca_components = max_pca_components
        self.targeted_variance_cv = targeted_variance_cv

    def build(self, depth_paths: list[str | Path]) -> PONModel:
        """Build a PON model from a list of depth-file paths."""

        samples = [
            DepthSample(sample=Path(path).stem, depth_file=Path(path)) for path in depth_paths
        ]
        return self.build_from_samples(samples)

    def build_from_samples(self, samples: list[DepthSample]) -> PONModel:
        """Build a PON model.

        Depths are divided by gene-specific tiling factors before comparison
        with the aggregate one-copy target baseline. Genes absent from the
        tiling file are treated as 1X.
        """

        if len(samples) < 2:
            raise ValueError("At least two normal depth files are required to build a PON")
        normalized_rows: list[pd.Series] = []
        pds_frames: list[pd.DataFrame] = []
        background_cvs: list[float] = []
        sample_stats: list[dict[str, object]] = []

        for sample in samples:
            profile = DepthProfile.from_path(sample.depth_file)
            region_means = profile.region_means(self.bed.regions)
            normalized = normalize_region_means(region_means, self.tiling_factors)
            normalized.name = sample.sample
            normalized_rows.append(normalized)
            adjusted_regions = add_region_tiling(region_means, self.tiling_factors)
            background_values = adjusted_regions.loc[
                adjusted_regions["tiling_factor"] == 1.0, "adjusted_depth"
            ].dropna()
            if not background_values.empty:
                mean_depth = float(background_values.mean())
                background_cv = float(background_values.std(ddof=1) / mean_depth)
                background_cvs.append(background_cv)
            else:
                mean_depth = float("nan")
                background_cv = float("nan")
            pds = profile.pds_signal(self.bed.pds_regions, self.bed.cyp_regions)
            pds_points = int(len(pds))
            if not pds.empty:
                pds = normalize_pds_signal(pds, region_means, self.bed.regions, self.tiling_factors)
                pds["normalized_depth"] = pds["depth"]
                pds["sample"] = sample.sample
                pds_frames.append(pds)
            sample_stats.append(
                {
                    "sample": sample.sample,
                    "depth_file": str(sample.depth_file),
                    "mean_1x_adjusted_depth": mean_depth,
                    "background_cv": background_cv,
                    "covered_regions": int(region_means["mean_depth"].notna().sum()),
                    "pds_points": pds_points,
                }
            )

        matrix = pd.DataFrame(normalized_rows)
        corrected_matrix, components = self._pca_correct(matrix)
        region_stats = self._region_stats(corrected_matrix)
        pds_stats = self._pds_stats(pds_frames)
        modality = self._infer_modality(background_cvs)
        pds_mode = "explicit_pds" if self.bed.pds_regions else "dense_cyp_target_fallback"

        return PONModel(
            version=__version__,
            bed_regions=self.bed.regions,
            region_stats=region_stats,
            pds_stats=pds_stats,
            modality=modality,
            pca_components=components,
            region_order=list(matrix.columns),
            pds_mode=pds_mode,
            tiling_factors=serializable_tiling_factors(self.tiling_factors),
            sample_stats=pd.DataFrame(sample_stats),
        )

    def _pca_correct(self, matrix: pd.DataFrame) -> tuple[pd.DataFrame, np.ndarray]:
        filled = matrix.copy()
        filled = filled.apply(lambda col: col.fillna(float(col.median())), axis=0)
        centered = filled - filled.mean(axis=0)
        if min(centered.shape) <= 1:
            return filled, np.empty((0, centered.shape[1]))
        _, singular_values, vt = np.linalg.svd(centered.to_numpy(dtype=float), full_matrices=False)
        n_components = min(self.max_pca_components, vt.shape[0] - 1, centered.shape[0] - 1)
        if n_components <= 0:
            return filled, np.empty((0, centered.shape[1]))
        components = vt[:n_components, :]
        scores = centered.to_numpy(dtype=float) @ components.T
        reconstruction = scores @ components
        corrected = filled - reconstruction
        corrected += filled.mean(axis=0)
        _ = singular_values
        return pd.DataFrame(corrected, index=matrix.index, columns=matrix.columns), components

    @staticmethod
    def _region_stats(matrix: pd.DataFrame) -> pd.DataFrame:
        pon_sd = matrix.std(axis=0, ddof=1).replace(0.0, np.nan).fillna(0.05)
        pon_sd = pon_sd.clip(lower=0.05)
        return pd.DataFrame(
            {
                "name": matrix.columns,
                "pon_mean": matrix.mean(axis=0).to_numpy(dtype=float),
                "pon_sd": pon_sd.to_numpy(dtype=float),
            }
        )

    @staticmethod
    def _pds_stats(pds_frames: list[pd.DataFrame]) -> pd.DataFrame:
        if not pds_frames:
            return pd.DataFrame(columns=["chrom", "pos0", "mean_depth", "sd_depth", "n"])
        pds_all = pd.concat(pds_frames, ignore_index=True)
        stats = (
            pds_all.groupby(["chrom", "pos0"])["normalized_depth"]
            .agg(["mean", "std", "count"])
            .reset_index()
            .rename(columns={"mean": "mean_depth", "std": "sd_depth", "count": "n"})
        )
        stats["sd_depth"] = stats["sd_depth"].replace(0.0, np.nan).fillna(0.05)
        stats["sd_depth"] = stats["sd_depth"].clip(lower=0.05)
        return stats

    def _infer_modality(self, cvs: list[float]) -> str:
        if not cvs:
            return "unknown"
        median_cv = float(np.nanmedian(cvs))
        return "wgs" if median_cv < self.targeted_variance_cv else "targeted_panel"
