"""CNV calling orchestration for ECHO."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd

from echo.algorithm import BPSDAlgorithm, ChangePoint
from echo.depth import DepthProfile
from echo.pon import PONModel
from echo.regions import BedRegion
from echo.tiling import normalize_pds_signal, normalize_region_means


@dataclass(frozen=True)
class CallArtifacts:
    """Intermediate and final results from one sample call."""

    report: dict[str, Any]
    cyp_region_calls: pd.DataFrame
    standard_region_calls: pd.DataFrame
    pds_signal: pd.DataFrame
    breakpoints: list[ChangePoint]


@dataclass(frozen=True)
class OutputPaths:
    """Files written for one sample call."""

    report: Path
    cnv_calls: Path
    gene_summary: Path
    cyp2d_plot: Path | None
    standard_gene_plots: dict[str, Path]


class CNVCaller:
    """Call CYP2D copy number and hybrid breakpoints from one sample."""

    def __init__(self, pon: PONModel, algorithm: BPSDAlgorithm | None = None) -> None:
        self.pon = pon
        self.algorithm = algorithm or BPSDAlgorithm()

    def call(
        self, depth_path: str | Path, sample: str | None = None, genes: list[str] | None = None
    ) -> dict[str, Any]:
        """Run ECHO on a single depth file.

        Parameters
        ----------
        depth_path
            Input per-base depth file.
        sample
            Optional sample identifier for the report.
        genes
            Optional requested gene symbols to call in addition to the CYP2D
            locus.

        Returns
        -------
        dict
            JSON-serializable report.
        """

        return self._call_sample(depth_path, sample, genes).report

    def write_outputs(
        self,
        depth_path: str | Path,
        output_dir: str | Path,
        sample: str,
        genes: list[str] | None = None,
        plot: bool = False,
        plot_format: str = "png",
    ) -> OutputPaths:
        """Call one sample and write requested outputs.

        Parameters
        ----------
        depth_path
            Input per-base depth file.
        sample
            Sample identifier used as the output filename prefix.
        output_dir
            Directory where ECHO writes the fixed output file set.
        genes
            Optional requested gene symbols to call and plot.
        plot
            Whether to write CYP2D and standard-gene plots.
        plot_format
            Plot file format: ``png``, ``svg``, or ``pdf``.

        Returns
        -------
        OutputPaths
            Paths written for the sample.
        """

        artifacts = self._call_sample(depth_path, sample, genes)
        report = artifacts.report
        standard_genes = [str(gene) for gene in report["standard_genes"]]
        output_directory = _resolve_output_dir(output_dir, sample)
        plot_extension = _normalize_plot_format(plot_format)
        paths = self._output_paths(
            output_directory, sample, standard_genes, plot=plot, plot_extension=plot_extension
        )
        paths.report.parent.mkdir(parents=True, exist_ok=True)

        self.write_report(report, paths.report)
        self.write_cnv_calls(report, paths.cnv_calls)
        self.write_gene_copy_numbers(report, paths.gene_summary)

        if not plot:
            self._remove_plot_outputs(output_directory, sample, standard_genes)
            return paths

        from echo.plotting import plot_cnv_call, plot_standard_gene_call

        if paths.cyp2d_plot is None:
            raise RuntimeError("CYP2D plot path is unavailable while plotting is enabled")
        plot_cnv_call(
            region_calls=artifacts.cyp_region_calls,
            pds_signal=artifacts.pds_signal,
            breakpoints=artifacts.breakpoints,
            output_path=paths.cyp2d_plot,
            sample=str(report["sample"]),
        )
        for gene, plot_path in paths.standard_gene_plots.items():
            gene_calls = artifacts.standard_region_calls[
                artifacts.standard_region_calls["gene"] == gene
            ]
            if gene_calls.empty:
                continue
            plot_standard_gene_call(
                region_calls=gene_calls,
                output_path=plot_path,
                sample=str(report["sample"]),
                gene=gene,
            )
        return paths

    def write_report(self, report: dict[str, Any], output_path: str | Path) -> None:
        """Write an ECHO report.

        JSON preserves the full structured report. TSV output uses the detailed
        CNV calls schema.
        """

        out = Path(output_path)
        if out.suffix.lower() == ".tsv":
            self.write_cnv_calls(report, out)
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

    def _call_sample(
        self,
        depth_path: str | Path,
        sample: str | None = None,
        genes: list[str] | None = None,
    ) -> CallArtifacts:
        """Call CYP2D and optional standard genes for one sample.

        Parameters
        ----------
        depth_path
            Input per-base depth file.
        sample
            Optional sample identifier.
        genes
            Optional standard gene symbols selected from the PON BED.

        Returns
        -------
        CallArtifacts
            Structured report plus plotting intermediates.
        """

        profile = DepthProfile.from_path(depth_path)
        bed = self.pon.bed
        region_means = profile.region_means(bed.regions)
        sample_norm = normalize_region_means(region_means, self.pon.tiling_factors)
        standard_genes = _normalize_gene_list(genes)
        selected_regions = self._selected_regions(standard_genes)
        region_calls = self._copy_number_by_region(sample_norm, selected_regions)
        cyp_names = {region.name for region in bed.cyp_regions}
        cyp_region_calls = region_calls[region_calls["name"].isin(cyp_names)].copy()
        standard_region_calls = region_calls[region_calls["gene"].isin(standard_genes)].copy()
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
            "standard_genes": standard_genes,
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
        return CallArtifacts(
            report=report,
            cyp_region_calls=cyp_region_calls,
            standard_region_calls=standard_region_calls,
            pds_signal=signal,
            breakpoints=breakpoints,
        )

    @staticmethod
    def _output_paths(
        output_dir: str | Path,
        sample: str,
        standard_genes: list[str],
        plot: bool,
        plot_extension: str,
    ) -> OutputPaths:
        """Return hardcoded ECHO output paths for one sample.

        Parameters
        ----------
        output_dir
            Output directory.
        sample
            Sample identifier used as filename prefix.
        standard_genes
            Optional standard genes requested for additional plotting.
        plot
            Whether plot paths are part of the output set.
        plot_extension
            Plot filename extension.

        Returns
        -------
        OutputPaths
            Fixed output paths.
        """

        out = Path(output_dir)
        standard_gene_plots = {
            gene: out / f"{sample}.{_filename_gene(gene)}.cnv.{plot_extension}"
            for gene in standard_genes
        }
        return OutputPaths(
            report=out / f"{sample}.echo.json",
            cnv_calls=out / f"{sample}.cnv.tsv",
            gene_summary=out / f"{sample}.genes.tsv",
            cyp2d_plot=out / f"{sample}.cyp2d6_cyp2d7.cnv.{plot_extension}" if plot else None,
            standard_gene_plots=standard_gene_plots if plot else {},
        )

    @staticmethod
    def _remove_plot_outputs(
        output_dir: str | Path, sample: str, standard_genes: list[str]
    ) -> None:
        """Remove fixed plot outputs when plotting is disabled.

        Parameters
        ----------
        output_dir
            Output directory.
        sample
            Sample identifier used as filename prefix.
        standard_genes
            Optional standard genes requested for additional plotting.

        Returns
        -------
        None
            Existing plot files are removed.
        """

        out = Path(output_dir)
        plot_paths = []
        for extension in ("png", "svg", "pdf"):
            plot_paths.append(out / f"{sample}.cyp2d6_cyp2d7.cnv.{extension}")
            plot_paths.extend(
                out / f"{sample}.{_filename_gene(gene)}.cnv.{extension}" for gene in standard_genes
            )
        for path in plot_paths:
            if path.exists():
                path.unlink()

    def _selected_regions(self, genes: list[str]) -> list[BedRegion]:
        """Return CYP2D and requested-gene regions.

        Parameters
        ----------
        genes
            Requested gene symbols.

        Returns
        -------
        list[BedRegion]
            Unique target regions used for copy-number calling.
        """

        bed = self.pon.bed
        selected: dict[str, BedRegion] = {region.name: region for region in bed.cyp_regions}
        if genes:
            matched = bed.regions_for_genes(genes)
            matched_genes = {region.gene for region in matched}
            missing = sorted(set(genes) - matched_genes)
            if missing:
                raise ValueError(
                    f"Requested genes are not present in the PON BED: {', '.join(missing)}"
                )
            selected.update({region.name: region for region in matched})
        return list(selected.values())

    def _copy_number_by_region(
        self, sample_norm: pd.Series, regions: list[BedRegion]
    ) -> pd.DataFrame:
        stats = self.pon.region_stats.set_index("name")
        rows = []
        for region in regions:
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
        columns = [
            "name",
            "gene",
            "feature",
            "exon",
            "chrom",
            "start",
            "end",
            "copy_number",
            "integer_copy_number",
            "z_score",
        ]
        return pd.DataFrame(rows).reindex(columns=columns)

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


def _normalize_gene_list(genes: list[str] | None) -> list[str]:
    """Return unique requested genes in input order.

    Parameters
    ----------
    genes
        Gene symbols from the CLI or API.

    Returns
    -------
    list[str]
        Unique, non-empty gene symbols.
    """

    normalized: list[str] = []
    seen: set[str] = set()
    for gene in genes or []:
        stripped = gene.strip()
        if not stripped or stripped in seen:
            continue
        normalized.append(stripped)
        seen.add(stripped)
    return normalized


def _filename_gene(gene: str) -> str:
    """Return a filesystem-safe gene token for output filenames.

    Parameters
    ----------
    gene
        Gene symbol.

    Returns
    -------
    str
        Lowercase filename token.
    """

    safe = "".join(char.lower() if char.isalnum() else "_" for char in gene.strip())
    return safe.strip("_") or "gene"


def _normalize_plot_format(plot_format: str) -> str:
    """Validate and normalize a plot format.

    Parameters
    ----------
    plot_format
        Requested plot format.

    Returns
    -------
    str
        Normalized extension.
    """

    normalized = plot_format.strip().lower()
    if normalized not in {"png", "svg", "pdf"}:
        raise ValueError("plot_format must be one of: png, svg, pdf")
    return normalized


def _resolve_output_dir(output_dir: str | Path, sample: str) -> Path:
    """Resolve output directory input.

    Parameters
    ----------
    output_dir
        Output directory path. If this points to a standard ECHO output file
        for the sample, the parent directory is used.
    sample
        Sample identifier used as filename prefix.

    Returns
    -------
    pathlib.Path
        Directory path for fixed ECHO outputs.
    """

    path = Path(output_dir)
    if _looks_like_echo_output_file(path, sample):
        return path.parent
    if path.exists() and path.is_file():
        raise ValueError(f"--output-dir must be a directory, not a file: {path}")
    return path


def _looks_like_echo_output_file(path: Path, sample: str) -> bool:
    """Return whether a path looks like a fixed ECHO output file.

    Parameters
    ----------
    path
        Candidate path.
    sample
        Sample identifier used as filename prefix.

    Returns
    -------
    bool
        True when the path is a known ECHO output filename.
    """

    if not path.name.startswith(f"{sample}."):
        return False
    known_suffixes = (
        ".echo.json",
        ".cnv.tsv",
        ".genes.tsv",
        ".cnv.png",
        ".cnv.svg",
        ".cnv.pdf",
    )
    return path.name.endswith(known_suffixes)
