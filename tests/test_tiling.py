from __future__ import annotations

from pathlib import Path

import pandas as pd

from echo.tiling import (
    add_region_tiling,
    dwbn_baseline_from_values,
    dwbn_stats_from_adjusted_regions,
    load_tiling_factors,
    normalize_region_means,
)


def test_tiling_normalization_uses_configured_gene_copy_factor(tmp_path: Path) -> None:
    """Configured 2X targets should normalize onto the same 1X baseline."""

    tiling_path = tmp_path / "tiling.tsv"
    tiling_path.write_text("gene\ttiling\nCYP2D6\t2\nCYP2D7\t2\n", encoding="utf-8")

    factors = load_tiling_factors(tiling_path)
    region_means = pd.DataFrame(
        [
            {"name": "A_ex1", "gene": "A", "mean_depth": 100.0},
            {"name": "CYP2D6_ex1", "gene": "CYP2D6", "mean_depth": 200.0},
            {"name": "CYP2D8P_ex1", "gene": "CYP2D8P", "mean_depth": 100.0},
        ]
    )

    normalized = normalize_region_means(region_means, factors)

    assert factors == {"CYP2D6": 2.0, "CYP2D7": 2.0}
    assert normalized["A_ex1"] == 1.0
    assert normalized["CYP2D6_ex1"] == 1.0
    assert normalized["CYP2D8P_ex1"] == 1.0


def test_dwbn_downweights_duplicated_background_region() -> None:
    region_means = pd.DataFrame(
        [
            {"name": "BG01", "gene": "A", "mean_depth": 100.0},
            {"name": "BG02", "gene": "B", "mean_depth": 101.0},
            {"name": "BG03", "gene": "C", "mean_depth": 99.0},
            {"name": "BG04", "gene": "D", "mean_depth": 100.0},
            {"name": "BG05", "gene": "E", "mean_depth": 100.0},
            {"name": "BG06", "gene": "F", "mean_depth": 98.0},
            {"name": "BG07", "gene": "G", "mean_depth": 102.0},
            {"name": "BG08", "gene": "H", "mean_depth": 100.0},
            {"name": "BG09", "gene": "I", "mean_depth": 101.0},
            {"name": "BG10", "gene": "J", "mean_depth": 99.0},
            {"name": "DUP", "gene": "K", "mean_depth": 200.0},
            {"name": "DUP2", "gene": "L", "mean_depth": 205.0},
        ]
    )

    adjusted = add_region_tiling(region_means, {})
    baseline = dwbn_stats_from_adjusted_regions(adjusted)
    normalized = normalize_region_means(region_means, {})

    assert 99.0 <= baseline.baseline <= 106.0
    assert baseline.raw_count == 12
    assert baseline.kde_mode < 110.0
    assert baseline.max_weight > baseline.min_weight
    assert 0.94 <= normalized["BG01"] <= 1.02
    assert normalized["DUP"] > 1.85


def test_dwbn_handles_single_background_region() -> None:
    baseline = dwbn_baseline_from_values(pd.Series([87.0]))

    assert baseline.baseline == 87.0
    assert baseline.kde_mode == 87.0
    assert baseline.raw_count == 1
