"""BED parsing and target-region utilities for ECHO."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

CYP_GENES = ("CYP2D6", "CYP2D7", "CYP2D8P")


@dataclass(frozen=True)
class BedRegion:
    """A single half-open BED interval.

    Parameters
    ----------
    chrom
        Chromosome name.
    start
        Zero-based inclusive start.
    end
        Zero-based exclusive end.
    name
        Annotation from the BED fourth column.
    gene
        Parsed gene symbol.
    feature
        Parsed feature token, such as ``ex9`` or ``in4``.
    """

    chrom: str
    start: int
    end: int
    name: str
    gene: str
    feature: str

    @property
    def length(self) -> int:
        """Return interval length in bases."""

        return self.end - self.start

    @property
    def is_cyp2d(self) -> bool:
        """Return whether the interval belongs to the CYP2D locus."""

        return any(self.gene.startswith(gene) for gene in CYP_GENES)

    @property
    def is_exon(self) -> bool:
        """Return whether the parsed feature is an exon."""

        return self.feature.startswith("ex") or "_ex" in self.name

    @property
    def exon_label(self) -> str:
        """Return a compact exon-level label."""

        match = re.search(r"(ex\d+[A-Za-z0-9]*)", self.name)
        return match.group(1) if match else self.feature


class AnnotatedBed:
    """Parsed annotated design BED.

    The parser derives gene names from the fourth column at runtime. No CYP2D
    coordinates are hardcoded.
    """

    def __init__(self, regions: list[BedRegion]) -> None:
        self.regions = regions

    @classmethod
    def from_path(cls, path: str | Path) -> AnnotatedBed:
        """Read an annotated BED file.

        Parameters
        ----------
        path
            BED path with at least four columns.

        Returns
        -------
        AnnotatedBed
            Parsed target design.
        """

        regions: list[BedRegion] = []
        with Path(path).open("r", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 4:
                    raise ValueError(f"BED line {line_number} has fewer than four columns")
                chrom, start_s, end_s, name = fields[:4]
                start = int(start_s)
                end = int(end_s)
                if end <= start:
                    raise ValueError(f"BED line {line_number} has non-positive interval length")
                gene, feature = parse_annotation(name)
                regions.append(BedRegion(chrom, start, end, name, gene, feature))
        if not regions:
            raise ValueError(f"No target regions found in {path}")
        return cls(regions)

    @property
    def cyp_regions(self) -> list[BedRegion]:
        """Return CYP2D6/CYP2D7/CYP2D8P intervals."""

        return [region for region in self.regions if region.is_cyp2d]

    @property
    def background_regions(self) -> list[BedRegion]:
        """Return non-CYP regions used for one-copy baseline estimation."""

        return [region for region in self.regions if not region.is_cyp2d]

    def regions_for_genes(self, genes: list[str]) -> list[BedRegion]:
        """Return regions matching requested gene symbols.

        Parameters
        ----------
        genes
            Gene symbols to select from the annotated BED.

        Returns
        -------
        list[BedRegion]
            Target regions whose parsed gene symbol matches one of the
            requested genes.
        """

        requested = {gene.strip() for gene in genes if gene.strip()}
        return [region for region in self.regions if region.gene in requested]

    @property
    def cyp_exons(self) -> list[BedRegion]:
        """Return CYP exon intervals."""

        return [region for region in self.cyp_regions if region.is_exon]

    @property
    def pds_regions(self) -> list[BedRegion]:
        """Return single-base CYP paralog-differentiating site intervals.

        The design may encode known PDS/SNP rows as 1 bp CYP2D intervals. When
        such rows are absent, callers can fall back to dense CYP target bases.
        """

        return [
            region
            for region in self.cyp_regions
            if region.length == 1
            and (
                "rs" in region.name.lower()
                or "pds" in region.name.lower()
                or "snp" in region.name.lower()
            )
        ]

    def to_frame(self) -> pd.DataFrame:
        """Return regions as a dataframe."""

        return pd.DataFrame(
            [
                {
                    "chrom": region.chrom,
                    "start": region.start,
                    "end": region.end,
                    "name": region.name,
                    "gene": region.gene,
                    "feature": region.feature,
                    "is_cyp2d": region.is_cyp2d,
                }
                for region in self.regions
            ]
        )


def parse_annotation(name: str) -> tuple[str, str]:
    """Parse a design BED annotation into gene and feature tokens."""

    for gene in CYP_GENES:
        if name.startswith(gene):
            feature = name.split("_")[-1] if "_" in name else name
            return gene, feature
    if "_" not in name:
        return name, name
    gene = name.split("_", maxsplit=1)[0]
    feature = name.split("_")[-1]
    return gene, feature
