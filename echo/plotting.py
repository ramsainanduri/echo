"""Scientific plotting for ECHO CNV calls."""

from __future__ import annotations

from pathlib import Path
from typing import Any, cast

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

from echo.algorithm import ChangePoint


def plot_cnv_call(
    region_calls: pd.DataFrame,
    pds_signal: pd.DataFrame,
    breakpoints: list[ChangePoint],
    output_path: str | Path,
    sample: str,
) -> None:
    """Create a CNV and BPSD diagnostic plot.

    The figure is intentionally denser than a simple CNV scatter plot. It shows
    absolute copy-number calls, PON-normalized z-scores, and the PDS depth-ratio
    signal used for Bayesian breakpoint detection.
    """

    out = Path(output_path)
    with plt.rc_context(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    ):
        fig, axes = plt.subplots(
            nrows=4,
            ncols=1,
            figsize=(14, 10),
            sharex=True,
            gridspec_kw={"height_ratios": [0.45, 1.05, 1.0, 1.35], "hspace": 0.18},
        )
        fig.suptitle(
            f"ECHO CYP2D CNV and hybrid profile: {sample}",
            fontsize=15,
            fontweight="bold",
            y=0.985,
        )
        fig.text(
            0.5,
            0.952,
            "Per-base depth normalized against PON; PDS shifts mark candidate CYP2D6/7 hybrids",
            fontsize=9,
            color="#4D4D4D",
            ha="center",
            va="top",
        )

        if not region_calls.empty:
            calls = _prepare_region_axis(region_calls)
            exon_ticks = _exon_ticks(calls)
            gene_spans = _gene_spans(calls)
            signal = _map_signal_to_region_axis(pds_signal, calls)
            mapped_breakpoints = _map_breakpoints_to_region_axis(breakpoints, calls)
            _plot_gene_model(axes[0], calls, gene_spans)
            _plot_copy_number_panel(axes[1], calls)
            z_scatter = _plot_z_score_panel(axes[2], calls)
            colorbar = fig.colorbar(z_scatter, ax=axes[2], pad=0.01, fraction=0.035)
            colorbar.set_label("z-score", rotation=270, labelpad=12)
        else:
            exon_ticks = ([], [])
            gene_spans = []
            signal = pds_signal
            mapped_breakpoints = [
                (breakpoint, float(index)) for index, breakpoint in enumerate(breakpoints)
            ]
            axes[0].text(0.5, 0.5, "No CYP regions evaluated", ha="center", va="center")
            axes[1].text(0.5, 0.5, "No copy-number calls available", ha="center", va="center")
            axes[2].text(0.5, 0.5, "No z-scores available", ha="center", va="center")

        _plot_pds_panel(axes[3], signal, mapped_breakpoints, exon_ticks)
        _apply_shared_style(axes)
        fig.align_ylabels(axes)
        fig.savefig(out, dpi=300, bbox_inches="tight", pad_inches=0.16)
        plt.close(fig)


def plot_standard_gene_call(
    region_calls: pd.DataFrame,
    output_path: str | Path,
    sample: str,
    gene: str,
) -> None:
    """Create a standard-gene CNV diagnostic plot.

    Parameters
    ----------
    region_calls
        Region-level copy-number calls for one requested gene. Required
        columns are ``feature``, ``exon``, ``start``, ``end``,
        ``copy_number``, and ``z_score``.
    output_path
        PNG path to write.
    sample
        Sample identifier shown in the figure title.
    gene
        Gene symbol shown in the figure title.

    Returns
    -------
    None
        The plot is written to ``output_path``.
    """

    out = Path(output_path)
    calls = _prepare_region_axis(region_calls)
    tick_positions, tick_labels = _standard_gene_ticks(calls)
    with plt.rc_context(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    ):
        fig, axes = plt.subplots(
            nrows=3,
            ncols=1,
            figsize=(10, 7),
            sharex=True,
            gridspec_kw={"height_ratios": [0.45, 1.05, 1.05], "hspace": 0.18},
        )
        fig.suptitle(
            f"ECHO CNV profile: {sample} {gene}",
            fontsize=14,
            fontweight="bold",
            y=0.985,
        )
        _plot_standard_gene_model(axes[0], calls, gene)
        _plot_standard_copy_number_panel(axes[1], calls, gene)
        z_scatter = _plot_z_score_panel(axes[2], calls)
        colorbar = fig.colorbar(z_scatter, ax=axes[2], pad=0.01, fraction=0.04)
        colorbar.set_label("z-score", rotation=270, labelpad=12)
        axes[2].set_xticks(tick_positions)
        axes[2].set_xticklabels(tick_labels, rotation=90, ha="center", va="top", fontsize=7)
        axes[2].set_xlabel("Target interval")
        _apply_shared_style(axes)
        fig.align_ylabels(axes)
        fig.savefig(out, dpi=300, bbox_inches="tight", pad_inches=0.16)
        plt.close(fig)


def _plot_standard_gene_model(axis: Any, calls: pd.DataFrame, gene: str) -> None:
    """Draw target blocks for one standard gene.

    Parameters
    ----------
    axis
        Matplotlib axis.
    calls
        Prepared region calls with ``plot_x``.
    gene
        Gene symbol.

    Returns
    -------
    None
        The axis is modified in place.
    """

    color = "#0072B2"
    axis.axhline(0.5, color="#666666", linewidth=0.8)
    for _, row in calls.iterrows():
        axis.add_patch(
            Rectangle(
                (float(row["plot_x"]) - 0.34, 0.28),
                0.68,
                0.44,
                facecolor=color,
                edgecolor="white",
                linewidth=0.6,
                alpha=0.9,
            )
        )
    axis.text(
        float(calls["plot_x"].median()),
        0.92,
        gene,
        ha="center",
        va="top",
        fontsize=9,
        fontweight="bold",
        color=color,
    )
    axis.set_ylim(0, 1)
    axis.set_yticks([])
    axis.set_ylabel("Gene\nmodel", rotation=0, ha="right", va="center")
    axis.set_xlim(-0.8, float(calls["plot_x"].max()) + 0.8)


def _plot_standard_copy_number_panel(axis: Any, calls: pd.DataFrame, gene: str) -> None:
    """Draw standard-gene copy number across target intervals.

    Parameters
    ----------
    axis
        Matplotlib axis.
    calls
        Prepared region calls with ``plot_x`` and ``copy_number``.
    gene
        Gene symbol used for the legend.

    Returns
    -------
    None
        The axis is modified in place.
    """

    color = "#0072B2"
    axis.plot(
        calls["plot_x"],
        calls["copy_number"],
        color=color,
        linestyle="--",
        linewidth=1.4,
        alpha=0.75,
        zorder=2,
    )
    axis.scatter(
        calls["plot_x"],
        calls["copy_number"],
        s=46,
        label=gene,
        color=color,
        edgecolor="white",
        linewidth=0.6,
        zorder=3,
    )
    for cn in (0, 1, 2, 3, 4):
        axis.axhline(cn, color="#D0D0D0", linewidth=0.7, zorder=0)
    axis.set_ylabel("Absolute CN")
    axis.set_ylim(-0.25, max(4.25, float(calls["copy_number"].max()) + 0.45))
    axis.set_yticks([0, 1, 2, 3, 4])
    axis.set_title("Copy-number state by target interval", loc="left", fontweight="bold", pad=8)
    axis.legend(loc="upper right", frameon=False)


def _plot_gene_model(
    axis: Any, calls: pd.DataFrame, gene_spans: list[tuple[str, float, float]]
) -> None:
    colors = _gene_colors()
    exons = calls[calls["feature"].astype(str).str.contains("ex", regex=False)].copy()
    axis.axhline(0.5, color="#666666", linewidth=0.8)
    for gene, group in exons.groupby("gene", sort=False):
        color = colors.get(str(gene), "#4D4D4D")
        for _, row in group.iterrows():
            axis.add_patch(
                Rectangle(
                    (float(row["plot_x"]) - 0.32, 0.28),
                    0.64,
                    0.44,
                    facecolor=color,
                    edgecolor="white",
                    linewidth=0.6,
                    alpha=0.9,
                )
            )
    for gene, start, end in gene_spans:
        axis.text(
            (start + end) / 2.0,
            0.92,
            gene,
            ha="center",
            va="top",
            fontsize=9,
            fontweight="bold",
            color=colors.get(gene, "#4D4D4D"),
        )
    axis.set_ylim(0, 1)
    axis.set_yticks([])
    axis.set_ylabel("Gene\nmodel", rotation=0, ha="right", va="center")
    axis.set_xlim(-0.8, float(calls["plot_x"].max()) + 0.8)


def _plot_copy_number_panel(axis: Any, calls: pd.DataFrame) -> None:
    colors = {"CYP2D6": "#0072B2", "CYP2D7": "#D55E00", "CYP2D8P": "#009E73"}
    for gene, group in calls.groupby("gene", sort=False):
        axis.scatter(
            group["plot_x"],
            group["copy_number"],
            s=42,
            label=str(gene),
            color=colors.get(str(gene), "#4D4D4D"),
            edgecolor="white",
            linewidth=0.5,
            zorder=3,
        )
        for _, row in group.iterrows():
            axis.hlines(
                y=float(row["copy_number"]),
                xmin=float(row["plot_x"]) - 0.42,
                xmax=float(row["plot_x"]) + 0.42,
                color=colors.get(str(gene), "#4D4D4D"),
                linewidth=2.4,
                alpha=0.75,
            )
    for cn in (0, 1, 2, 3, 4):
        axis.axhline(cn, color="#D0D0D0", linewidth=0.7, zorder=0)
    axis.set_ylabel("Absolute CN")
    axis.set_ylim(-0.25, max(4.25, float(calls["copy_number"].max()) + 0.45))
    axis.set_yticks([0, 1, 2, 3, 4])
    axis.set_title("Copy-number state by target interval", loc="left", fontweight="bold", pad=8)
    axis.legend(loc="upper right", frameon=False, ncols=3)


def _plot_z_score_panel(axis: Any, calls: pd.DataFrame) -> Any:
    scatter = axis.scatter(
        calls["plot_x"],
        calls["z_score"],
        s=40,
        c=calls["z_score"],
        cmap="coolwarm",
        vmin=-4,
        vmax=4,
        edgecolor="#4D4D4D",
        linewidth=0.35,
        zorder=3,
    )
    axis.axhline(0, color="#333333", linewidth=0.9)
    for threshold in (-3, -2, 2, 3):
        axis.axhline(
            threshold,
            color="#A0A0A0",
            linestyle="--" if abs(threshold) == 2 else ":",
            linewidth=0.8,
        )
    axis.set_ylabel("PON z-score")
    axis.set_title("Deviation from panel of normals", loc="left", fontweight="bold", pad=8)
    axis.text(
        0.01,
        0.94,
        "z = (sample depth - PON mean) / PON SD; negative values indicate lower depth",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=7.5,
        color="#4D4D4D",
        bbox={
            "boxstyle": "round,pad=0.25",
            "fc": "white",
            "ec": "#BDBDBD",
            "lw": 0.5,
            "alpha": 0.88,
        },
    )
    limit = max(4.0, min(12.0, float(np.nanmax(np.abs(calls["z_score"]))) + 0.75))
    axis.set_ylim(-limit, limit)
    return scatter


def _plot_pds_panel(
    axis: Any,
    pds_signal: pd.DataFrame,
    breakpoints: list[tuple[ChangePoint, float]],
    exon_ticks: tuple[list[float], list[str]],
) -> None:
    if pds_signal.empty:
        axis.text(0.5, 0.5, "No PDS signal available", ha="center", va="center")
        axis.set_ylabel("PDS ratio")
        return

    signal = pds_signal.sort_values("plot_x") if "plot_x" in pds_signal.columns else pds_signal
    axis.scatter(
        signal["plot_x"],
        signal["pds_ratio"],
        s=34,
        c=signal["z_score"] if "z_score" in signal.columns else signal["pds_ratio"],
        cmap="coolwarm",
        vmin=-4,
        vmax=4,
        alpha=0.95,
        edgecolor="#333333",
        linewidth=0.3,
        label="PDS points",
        zorder=4,
    )
    rolling_window = min(15, max(5, len(signal) // 10))
    rolling = (
        signal["pds_ratio"].rolling(window=rolling_window, center=True, min_periods=3).median()
    )
    axis.plot(
        signal["plot_x"], rolling, color="#E69F00", linewidth=2.4, label="Rolling median", zorder=5
    )
    axis.axhline(1.0, color="#333333", linewidth=0.9)
    for ratio in (0.5, 1.5):
        axis.axhline(ratio, color="#BDBDBD", linestyle="--", linewidth=0.8)
    if np.isfinite(signal["pds_ratio"]).any():
        lower = max(0.0, float(np.nanpercentile(signal["pds_ratio"], 2)) - 0.15)
        upper = max(float(np.nanpercentile(signal["pds_ratio"], 98)) + 0.15, 1.65)
        axis.set_ylim(lower, upper)
    else:
        upper = 1.65
    strongest_breakpoints = sorted(
        breakpoints, key=lambda item: item[0].log_bayes_factor, reverse=True
    )[:4]
    strongest_breakpoints = sorted(strongest_breakpoints, key=lambda item: item[1])
    label_offsets = [(5, -8), (5, -30), (5, -52), (5, -74)]
    label_y = upper - 0.03
    for index, (breakpoint, plot_x) in enumerate(strongest_breakpoints, start=1):
        axis.axvline(plot_x, color="#CC0000", linewidth=1.3, alpha=0.82, zorder=2)
        axis.annotate(
            f"BP{index}\nBF {breakpoint.log_bayes_factor:.0f}",
            xy=(plot_x, label_y),
            xytext=label_offsets[(index - 1) % len(label_offsets)],
            textcoords="offset points",
            va="top",
            ha="left",
            fontsize=7.5,
            color="#CC0000",
            bbox={
                "boxstyle": "round,pad=0.2",
                "fc": "white",
                "ec": "#CC0000",
                "lw": 0.6,
                "alpha": 0.9,
            },
        )
    axis.set_ylabel("PDS ratio")
    axis.set_title(
        "Paralog signal deconvolution and hybrid breakpoint support",
        loc="left",
        fontweight="bold",
        pad=8,
    )
    axis.text(
        0.01,
        0.97,
        "Red lines = Bayesian change points; BP = breakpoint candidate; BF = Bayes factor support",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=7.5,
        color="#4D4D4D",
        bbox={
            "boxstyle": "round,pad=0.25",
            "fc": "white",
            "ec": "#BDBDBD",
            "lw": 0.5,
            "alpha": 0.9,
        },
    )
    tick_positions, tick_labels = exon_ticks
    if tick_positions:
        axis.set_xticks(tick_positions)
        axis.set_xticklabels(tick_labels, rotation=90, ha="center", va="top", fontsize=7)
    axis.set_xlabel("CYP2D exon position")
    axis.legend(loc="upper right", frameon=False)


def _exon_ticks(calls: pd.DataFrame) -> tuple[list[float], list[str]]:
    exons = calls[calls["feature"].astype(str).str.contains("ex", regex=False)].copy()
    if exons.empty:
        return [], []
    grouped = (
        exons.groupby(["gene", "exon"], as_index=False)
        .agg(plot_x=("plot_x", "median"), start=("start", "min"))
        .sort_values("start")
    )
    positions = [float(value) for value in grouped["plot_x"].tolist()]
    labels = [
        f"{str(row.gene).replace('CYP', '')} {row.exon}" for row in grouped.itertuples(index=False)
    ]
    return positions, labels


def _standard_gene_ticks(calls: pd.DataFrame) -> tuple[list[float], list[str]]:
    """Return x-axis ticks for standard-gene plots.

    Parameters
    ----------
    calls
        Prepared region calls with ``plot_x``, ``exon``, ``feature``, and
        ``start`` columns.

    Returns
    -------
    tuple[list[float], list[str]]
        Tick positions and labels.
    """

    positions = [float(value) for value in calls["plot_x"].tolist()]
    labels: list[str] = []
    for index, row in enumerate(calls.itertuples(index=False), start=1):
        exon = str(getattr(row, "exon", "") or getattr(row, "feature", ""))
        labels.append(exon if exon and exon != "nan" else f"target{index}")
    return positions, labels


def _prepare_region_axis(region_calls: pd.DataFrame) -> pd.DataFrame:
    calls = region_calls.sort_values(["chrom", "start", "end"]).reset_index(drop=True).copy()
    calls["plot_x"] = np.arange(len(calls), dtype=float)
    return calls


def _map_signal_to_region_axis(pds_signal: pd.DataFrame, calls: pd.DataFrame) -> pd.DataFrame:
    if pds_signal.empty:
        return pds_signal
    mapped = pds_signal.copy()
    mapped["plot_x"] = [
        _coordinate_to_plot_x(int(cast(Any, row.pos0)), calls)
        for row in mapped.itertuples(index=False)
    ]
    return mapped.dropna(subset=["plot_x"])


def _map_breakpoints_to_region_axis(
    breakpoints: list[ChangePoint], calls: pd.DataFrame
) -> list[tuple[ChangePoint, float]]:
    mapped: list[tuple[ChangePoint, float]] = []
    for breakpoint in breakpoints:
        plot_x = _coordinate_to_plot_x(breakpoint.coordinate, calls)
        if np.isfinite(plot_x):
            mapped.append((breakpoint, float(plot_x)))
    return mapped


def _coordinate_to_plot_x(coordinate: int, calls: pd.DataFrame) -> float:
    containing = calls[(calls["start"] <= coordinate) & (calls["end"] > coordinate)]
    if not containing.empty:
        row = containing.iloc[0]
        width = max(float(row["end"] - row["start"]), 1.0)
        fraction = min(max((coordinate - float(row["start"])) / width, 0.0), 1.0)
        return float(row["plot_x"]) - 0.42 + 0.84 * fraction
    starts = calls["start"].to_numpy(dtype=float)
    plot_x = calls["plot_x"].to_numpy(dtype=float)
    return float(np.interp(coordinate, starts, plot_x))


def _gene_colors() -> dict[str, str]:
    return {"CYP2D6": "#0072B2", "CYP2D7": "#D55E00", "CYP2D8P": "#009E73"}


def _gene_spans(calls: pd.DataFrame) -> list[tuple[str, float, float]]:
    spans: list[tuple[str, float, float]] = []
    for gene, group in calls.groupby("gene", sort=False):
        spans.append(
            (str(gene), float(group["plot_x"].min()) - 0.5, float(group["plot_x"].max()) + 0.5)
        )
    return spans


def _apply_shared_style(axes: np.ndarray) -> None:
    for axis in axes:
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.spines["left"].set_color("#808080")
        axis.spines["bottom"].set_color("#808080")
        axis.grid(axis="y", color="#E6E6E6", linewidth=0.7)
        axis.tick_params(axis="both", colors="#333333", length=3)
