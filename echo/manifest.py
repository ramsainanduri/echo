"""Sample manifest parsing for ECHO."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass(frozen=True)
class DepthSample:
    """A sample identifier paired with a per-base depth file."""

    sample: str
    depth_file: Path


def load_depth_manifest(path: str | Path) -> list[DepthSample]:
    """Load a sample/depth-file manifest.

    Parameters
    ----------
    path
        CSV or TSV with required columns ``sample`` and ``depth_file``.
        Relative depth-file paths are resolved relative to the manifest file.

    Returns
    -------
    list of DepthSample
        Ordered manifest records.
    """

    manifest_path = Path(path)
    frame = pd.read_csv(manifest_path, sep=None, engine="python", comment="#")
    frame.columns = [str(column).strip() for column in frame.columns]
    required = {"sample", "depth_file"}
    missing = required - set(frame.columns)
    if missing:
        missing_s = ", ".join(sorted(missing))
        raise ValueError(f"Depth manifest is missing required column(s): {missing_s}")

    samples: list[DepthSample] = []
    seen: set[str] = set()
    for row_offset, (_, row) in enumerate(frame.iterrows(), start=2):
        sample = str(row["sample"]).strip()
        depth_value = str(row["depth_file"]).strip()
        line_number = row_offset
        if not sample:
            raise ValueError(f"Depth manifest line {line_number} has an empty sample")
        if sample in seen:
            raise ValueError(f"Depth manifest contains duplicate sample: {sample}")
        if not depth_value:
            raise ValueError(f"Depth manifest line {line_number} has an empty depth_file")
        depth_file = Path(depth_value)
        if not depth_file.is_absolute():
            depth_file = manifest_path.parent / depth_file
        if not depth_file.exists():
            raise FileNotFoundError(
                f"Depth file for sample {sample!r} does not exist: {depth_file}"
            )
        seen.add(sample)
        samples.append(DepthSample(sample=sample, depth_file=depth_file))

    if not samples:
        raise ValueError("Depth manifest contains no samples")
    return samples
