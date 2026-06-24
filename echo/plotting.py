"""Scientific plotting for ECHO CNV calls."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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
    fig, axes = plt.subplots(
        nrows=3,
        ncols=1,
        figsize=(13, 9),
        sharex=True,
        gridspec_kw={"height_ratios": [1.1, 1.0, 1.2]},
    )
    fig.suptitle(f"ECHO CNV profile: {sample}", fontsize=14, fontweight="bold")

    if not region_calls.empty:
        calls = region_calls.sort_values(["chrom", "start"]).copy()
        calls["midpoint"] = (calls["start"] + calls["end"]) / 2.0
        _plot_copy_number_panel(axes[0], calls)
        _plot_z_score_panel(axes[1], calls)
    else:
        axes[0].text(0.5, 0.5, "No CYP regions evaluated", ha="center", va="center")
        axes[1].text(0.5, 0.5, "No z-scores available", ha="center", va="center")

    _plot_pds_panel(axes[2], pds_signal, breakpoints)
    _apply_shared_style(axes)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out, dpi=220)
    plt.close(fig)


def _plot_copy_number_panel(axis: Any, calls: pd.DataFrame) -> None:
    colors = {"CYP2D6": "#0072B2", "CYP2D7": "#D55E00", "CYP2D8P": "#009E73"}
    for gene, group in calls.groupby("gene", sort=False):
        axis.scatter(
            group["midpoint"],
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
                xmin=float(row["start"]),
                xmax=float(row["end"]),
                color=colors.get(str(gene), "#4D4D4D"),
                linewidth=2.4,
                alpha=0.75,
            )
    for cn in (0, 1, 2, 3, 4):
        axis.axhline(cn, color="#D0D0D0", linewidth=0.7, zorder=0)
    axis.set_ylabel("Absolute CN")
    axis.set_ylim(bottom=-0.25)
    axis.legend(loc="upper right", frameon=False, ncols=3)


def _plot_z_score_panel(axis: Any, calls: pd.DataFrame) -> None:
    axis.scatter(
        calls["midpoint"],
        calls["z_score"],
        s=34,
        c=calls["z_score"],
        cmap="coolwarm",
        vmin=-4,
        vmax=4,
        edgecolor="white",
        linewidth=0.4,
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


def _plot_pds_panel(axis: Any, pds_signal: pd.DataFrame, breakpoints: list[ChangePoint]) -> None:
    if pds_signal.empty:
        axis.text(0.5, 0.5, "No PDS signal available", ha="center", va="center")
        axis.set_ylabel("PDS ratio")
        return

    signal = pds_signal.sort_values(["chrom", "pos0"])
    axis.plot(
        signal["pos0"],
        signal["pds_ratio"],
        color="#222222",
        linewidth=1.0,
        alpha=0.8,
        label="PDS ratio",
    )
    axis.scatter(
        signal["pos0"],
        signal["pds_ratio"],
        s=8,
        color="#6A51A3",
        alpha=0.45,
        zorder=3,
    )
    rolling = signal["pds_ratio"].rolling(window=25, center=True, min_periods=5).median()
    axis.plot(signal["pos0"], rolling, color="#E69F00", linewidth=2.0, label="Rolling median")
    axis.axhline(1.0, color="#333333", linewidth=0.9)
    for breakpoint in breakpoints:
        axis.axvline(breakpoint.coordinate, color="#CC0000", linewidth=1.4, alpha=0.85)
        axis.text(
            breakpoint.coordinate,
            float(np.nanmax(signal["pds_ratio"])),
            f"BF={breakpoint.log_bayes_factor:.1f}",
            rotation=90,
            va="top",
            ha="right",
            fontsize=8,
            color="#CC0000",
        )
    axis.set_ylabel("PDS ratio")
    axis.set_xlabel("Genomic coordinate (0-based)")
    axis.legend(loc="upper right", frameon=False)


def _apply_shared_style(axes: np.ndarray) -> None:
    for axis in axes:
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.grid(axis="y", color="#E6E6E6", linewidth=0.7)
