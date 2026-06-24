"""CNV calling orchestration for ECHO."""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd

from echo.algorithm import BPSDAlgorithm, ChangePoint
from echo.depth import DepthProfile
from echo.pon import PONModel
from echo.tiling import normalize_pds_signal, normalize_region_means


class CNVCaller:
    """Call CYP2D copy number and hybrid breakpoints from one sample."""

    def __init__(self, pon: PONModel, algorithm: BPSDAlgorithm | None = None) -> None:
        self.pon = pon
        self.algorithm = algorithm or BPSDAlgorithm()

    def call(self, depth_path: str | Path, sample: str | None = None) -> dict[str, Any]:
        """Run ECHO on a single depth file.

        Parameters
        ----------
        depth_path
            Input per-base depth file.

        Returns
        -------
        dict
            JSON-serializable report.
        """

        report, _, _, _ = self._call_with_intermediates(depth_path, sample)
        return report

    def write_outputs(
        self,
        depth_path: str | Path,
        report_path: str | Path,
        cnv_output_path: str | Path | None = None,
        gene_output_path: str | Path | None = None,
        plot_path: str | Path | None = None,
        sample: str | None = None,
    ) -> None:
        """Call one sample and write requested report, CNV, gene, and plot outputs."""

        report, region_calls, signal, breakpoints = self._call_with_intermediates(
            depth_path, sample
        )
        self.write_report(report, report_path)
        if cnv_output_path is not None:
            self.write_cnv_calls(report, cnv_output_path)
        if gene_output_path is not None:
            self.write_gene_copy_numbers(report, gene_output_path)
        if plot_path is not None:
            from echo.plotting import plot_cnv_call

            plot_cnv_call(
                region_calls=region_calls,
                pds_signal=signal,
                breakpoints=breakpoints,
                output_path=plot_path,
                sample=str(report["sample"]),
            )

    def write_report(self, report: dict[str, Any], output_path: str | Path) -> None:
        """Write a JSON or TSV ECHO report."""

        out = Path(output_path)
        if out.suffix.lower() == ".tsv":
            self._write_tsv(report, out)
            return
        with out.open("w", encoding="utf-8") as handle:
            json.dump(report, handle, indent=2, sort_keys=True)
            handle.write("\n")

    def write_cnv_calls(self, report: dict[str, Any], output_path: str | Path) -> None:
        """Write a dedicated CNV calls table."""

        rows: list[dict[str, Any]] = []
        for gene, call in report["gene_copy_number"].items():
            rows.append(
                {
                    "call_type": "gene",
                    "gene": gene,
                    "copy_number": call["copy_number"],
                    "integer_copy_number": call["integer_copy_number"],
                    "hybrid": "no",
                    "segments": call["segments"],
                }
            )
        for row in report["exon_copy_number"]:
            rows.append(
                {
                    "call_type": "exon",
                    "gene": row["gene"],
                    "copy_number": row["copy_number"],
                    "integer_copy_number": row["integer_copy_number"],
                    "hybrid": "no",
                    "feature": row["exon"],
                    "chrom": row["chrom"],
                    "start": row["start"],
                    "end": row["end"],
                    "z_score": row["z_score"],
                }
            )
        for row in report["hybrid_calls"]:
            rows.append(
                {
                    "call_type": "hybrid",
                    "gene": "CYP2D6/CYP2D7",
                    "hybrid": "yes",
                    "feature": row["feature"],
                    "chrom": row["chrom"],
                    "breakpoint_coordinate": row["coordinate"],
                    "log_bayes_factor": row["log_bayes_factor"],
                    "left_pds_ratio": row["left_pds_ratio"],
                    "right_pds_ratio": row["right_pds_ratio"],
                }
            )
        columns = [
            "call_type",
            "gene",
            "copy_number",
            "integer_copy_number",
            "hybrid",
            "feature",
            "chrom",
            "start",
            "end",
            "z_score",
            "breakpoint_coordinate",
            "log_bayes_factor",
            "left_pds_ratio",
            "right_pds_ratio",
            "segments",
        ]
        pd.DataFrame(rows).reindex(columns=columns).to_csv(output_path, sep="\t", index=False)

    def write_gene_copy_numbers(self, report: dict[str, Any], output_path: str | Path) -> None:
        """Write a compact gene-level copy-number table."""

        rows: list[dict[str, Any]] = []
        for gene, call in report["gene_copy_number"].items():
            copy_number = float(call["copy_number"])
            integer_copy_number = int(call["integer_copy_number"])
            rows.append(
                {
                    "gene": gene,
                    "CN": integer_copy_number,
                    "copy_number": copy_number,
                }
            )
        pd.DataFrame(rows).reindex(columns=["gene", "CN", "copy_number"]).to_csv(
            output_path, sep="\t", index=False
        )

    def _call_with_intermediates(
        self, depth_path: str | Path, sample: str | None = None
    ) -> tuple[dict[str, Any], pd.DataFrame, pd.DataFrame, list[ChangePoint]]:
        profile = DepthProfile.from_path(depth_path)
        bed = self.pon.bed
        region_means = profile.region_means(bed.regions)
        sample_norm = normalize_region_means(region_means, self.pon.tiling_factors)
        region_calls = self._copy_number_by_region(sample_norm)
        exon_calls = self._exon_calls(region_calls)
        gene_calls = self._gene_calls(region_calls)

        pds = profile.pds_signal(bed.pds_regions, bed.cyp_regions)
        pds = normalize_pds_signal(pds, region_means, bed.regions, self.pon.tiling_factors)
        signal, breakpoints = self.algorithm.deconvolve(pds, self.pon.pds_stats)
        breakpoint_rows = [
            {
                "chrom": str(signal.iloc[min(cp.right_index, len(signal) - 1)]["chrom"])
                if not signal.empty
                else None,
                **asdict(cp),
            }
            for cp in breakpoints
        ]
        hybrid_calls = self._hybrid_calls(breakpoint_rows)
        sample_id = sample or Path(depth_path).stem
        report = {
            "sample": sample_id,
            "pon_version": self.pon.version,
            "sequencing_modality": self.pon.modality,
            "pds_mode": self.pon.pds_mode,
            "tiling_factors": self.pon.tiling_factors,
            "gene_copy_number": gene_calls,
            "exon_copy_number": exon_calls,
            "breakpoints": breakpoint_rows,
            "hybrid_calls": hybrid_calls,
            "quality": {
                "regions_evaluated": int(len(region_calls)),
                "pds_points_evaluated": int(len(signal)),
                "median_abs_pds_z": float(np.nanmedian(np.abs(signal["z_score"])))
                if not signal.empty
                else None,
            },
        }
        return report, region_calls, signal, breakpoints

    def _copy_number_by_region(self, sample_norm: pd.Series) -> pd.DataFrame:
        stats = self.pon.region_stats.set_index("name")
        rows = []
        for region in self.pon.bed.cyp_regions:
            if region.name not in sample_norm.index or region.name not in stats.index:
                continue
            observed = float(sample_norm.loc[region.name])
            pon_mean = float(cast(Any, stats.loc[region.name, "pon_mean"]))
            pon_sd = float(cast(Any, stats.loc[region.name, "pon_sd"]))
            ratio = observed / max(pon_mean, 1e-12)
            copy_number = max(0.0, 2.0 * ratio)
            z_score = (observed - pon_mean) / max(pon_sd, 1e-12)
            rows.append(
                {
                    "name": region.name,
                    "gene": region.gene,
                    "feature": region.feature,
                    "exon": region.exon_label,
                    "chrom": region.chrom,
                    "start": region.start,
                    "end": region.end,
                    "copy_number": copy_number,
                    "integer_copy_number": int(np.rint(copy_number)),
                    "z_score": float(z_score),
                }
            )
        return pd.DataFrame(rows)

    @staticmethod
    def _exon_calls(region_calls: pd.DataFrame) -> list[dict[str, Any]]:
        if region_calls.empty:
            return []
        exons = region_calls[region_calls["feature"].str.contains("ex", regex=False)].copy()
        if exons.empty:
            return []
        grouped = (
            exons.groupby(["gene", "exon", "chrom", "start", "end"], as_index=False)
            .agg(copy_number=("copy_number", "median"), z_score=("z_score", "median"))
            .sort_values(["gene", "start"])
        )
        grouped["integer_copy_number"] = np.rint(grouped["copy_number"]).astype(int)
        return cast(list[dict[str, Any]], grouped.to_dict(orient="records"))

    @staticmethod
    def _gene_calls(region_calls: pd.DataFrame) -> dict[str, Any]:
        if region_calls.empty:
            return {}
        calls: dict[str, Any] = {}
        for gene, group in region_calls.groupby("gene"):
            cn = float(np.nanmedian(group["copy_number"]))
            calls[str(gene)] = {
                "copy_number": cn,
                "integer_copy_number": int(np.rint(cn)),
                "segments": int(len(group)),
            }
        return calls

    def _hybrid_calls(self, breakpoints: list[dict[str, Any]]) -> list[dict[str, Any]]:
        calls: list[dict[str, Any]] = []
        for row in breakpoints:
            coordinate = row.get("coordinate")
            if coordinate is None:
                continue
            region = self._region_at_coordinate(int(coordinate))
            calls.append(
                {
                    "event": "CYP2D6/CYP2D7 hybrid candidate",
                    "chrom": row.get("chrom"),
                    "coordinate": coordinate,
                    "feature": region.name if region is not None else None,
                    "log_bayes_factor": row.get("log_bayes_factor"),
                    "left_pds_ratio": row.get("left_mean"),
                    "right_pds_ratio": row.get("right_mean"),
                }
            )
        return calls

    def _region_at_coordinate(self, coordinate: int) -> Any | None:
        for region in self.pon.bed.cyp_regions:
            if region.start <= coordinate < region.end:
                return region
        return None

    def _write_tsv(self, report: dict[str, Any], output_path: Path) -> None:
        rows: list[dict[str, Any]] = []
        for row in report["exon_copy_number"]:
            rows.append({"record_type": "exon", **row})
        for row in report["breakpoints"]:
            rows.append({"record_type": "breakpoint", **row})
        pd.DataFrame(rows).to_csv(output_path, sep="\t", index=False)
