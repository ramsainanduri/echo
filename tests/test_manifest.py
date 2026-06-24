from __future__ import annotations

from pathlib import Path

import pytest

from echo.manifest import load_depth_manifest


def test_pon_manifest_loads_sample_depth_inputs(tmp_path: Path) -> None:
    """PON manifests should preserve sample IDs and resolve relative depth paths."""

    depth_dir = tmp_path / "depths"
    depth_dir.mkdir()
    depth_path = depth_dir / "HG001.depth.tsv"
    depth_path.write_text("chr1\t1\t30\n", encoding="utf-8")
    manifest = tmp_path / "pon_samples.tsv"
    manifest.write_text(
        "sample\tdepth_file\nHG001\tdepths/HG001.depth.tsv\n",
        encoding="utf-8",
    )

    samples = load_depth_manifest(manifest)

    assert samples[0].sample == "HG001"
    assert samples[0].depth_file == depth_path


def test_pon_manifest_rejects_duplicate_sample_ids(tmp_path: Path) -> None:
    """Duplicate normal IDs would make PON matrices ambiguous."""

    depth_path = tmp_path / "HG001.depth.tsv"
    depth_path.write_text("chr1\t1\t30\n", encoding="utf-8")
    manifest = tmp_path / "pon_samples.csv"
    manifest.write_text(
        "sample,depth_file\nHG001,HG001.depth.tsv\nHG001,HG001.depth.tsv\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="duplicate sample"):
        load_depth_manifest(manifest)
