from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from echo.algorithm import ChangePoint
from echo.plotting import plot_cnv_call
from echo.pon import PONModel
from echo.regions import BedRegion


def test_cnv_diagnostic_plot_is_written(tmp_path: Path) -> None:
    region_calls = pd.DataFrame(
        [
            {
                "chrom": "chr22",
                "gene": "CYP2D6",
                "start": 42126573,
                "end": 42126752,
                "copy_number": 2.0,
                "z_score": 0.1,
            },
            {
                "chrom": "chr22",
                "gene": "CYP2D7",
                "start": 42140203,
                "end": 42140456,
                "copy_number": 1.0,
                "z_score": -2.4,
            },
        ]
    )
    pds_signal = pd.DataFrame(
        {
            "chrom": ["chr22"] * 30,
            "pos0": np.arange(42128780, 42128810),
            "pds_ratio": [1.0] * 15 + [1.45] * 15,
            "z_score": [0.0] * 15 + [3.0] * 15,
        }
    )
    breakpoints = [
        ChangePoint(
            left_index=14,
            right_index=15,
            coordinate=42128795,
            log_bayes_factor=12.0,
            left_mean=1.0,
            right_mean=1.45,
        )
    ]
    output = tmp_path / "HG001.cnv.png"

    plot_cnv_call(region_calls, pds_signal, breakpoints, output, sample="HG001")

    assert output.exists()
    assert output.stat().st_size > 0


def test_pon_stats_output_contains_build_qc(tmp_path: Path) -> None:
    model = PONModel(
        version="0.1.0",
        bed_regions=[BedRegion("chr22", 1, 2, "CYP2D6_NM_000106.6_ex1", "CYP2D6", "ex1")],
        region_stats=pd.DataFrame(
            {"name": ["CYP2D6_NM_000106.6_ex1"], "pon_mean": [1.0], "pon_sd": [0.1]}
        ),
        pds_stats=pd.DataFrame(
            {"chrom": ["chr22"], "pos0": [1], "mean_depth": [1.0], "sd_depth": [0.1], "n": [2]}
        ),
        modality="targeted_panel",
        pca_components=np.empty((0, 1)),
        region_order=["CYP2D6_NM_000106.6_ex1"],
        pds_mode="explicit_pds",
        tiling_factors={"CYP2D6": 2.0},
        sample_stats=pd.DataFrame(
            {
                "sample": ["HG001"],
                "mean_1x_adjusted_depth": [100.0],
                "background_cv": [0.12],
                "covered_regions": [1],
                "pds_points": [1],
            }
        ),
    )
    output = tmp_path / "pon.stats.tsv"

    model.write_stats(output)

    stats = pd.read_csv(output, sep="\t")
    assert stats.loc[0, "sample"] == "HG001"
    assert stats.loc[0, "background_cv"] == 0.12
