from __future__ import annotations

from pathlib import Path

from echo.regions import AnnotatedBed


def test_annotated_target_design_classifies_cyp_locus_regions(tmp_path: Path) -> None:
    """The target design BED is the reference for gene and feature ownership."""

    bed_path = tmp_path / "target_design.bed"
    bed_path.write_text(
        "chr1\t0\t10\tGENE_NM_1_ex1\n"
        "chr22\t20\t30\tCYP2D6_NM_000106.6_ex9\n"
        "chr22\t40\t41\tCYP2D7_NR_002570.6_rsX\n",
        encoding="utf-8",
    )

    bed = AnnotatedBed.from_path(bed_path)

    assert len(bed.background_regions) == 1
    assert len(bed.cyp_regions) == 2
    assert len(bed.pds_regions) == 1
