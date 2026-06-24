from __future__ import annotations

from pathlib import Path

import pandas as pd

from echo.tiling import load_tiling_factors, normalize_region_means


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
