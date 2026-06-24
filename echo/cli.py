"""Command-line interface for ECHO."""

from __future__ import annotations

import argparse
from pathlib import Path

from echo.caller import CNVCaller
from echo.manifest import DepthSample, load_depth_manifest
from echo.pon import PanelOfNormalsBuilder, PONModel


def build_parser() -> argparse.ArgumentParser:
    """Create the CLI argument parser."""

    parser = argparse.ArgumentParser(
        prog="echo",
        description="ECHO: Bayesian paralog signal deconvolution for CYP2D CNV calling.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    build = subparsers.add_parser("build-pon", help="Build a panel-of-normals model.")
    build.add_argument(
        "--manifest",
        type=Path,
        help="CSV/TSV with columns sample and depth_file for PON normal samples.",
    )
    build.add_argument(
        "--depth-dir",
        type=Path,
        help="Directory of normal depth files. Used only when --manifest is not provided.",
    )
    build.add_argument("--bed", required=True, type=Path, help="Annotated design BED.")
    build.add_argument(
        "--tiling",
        type=Path,
        help=("Optional two-column gene tiling file. Genes not listed are treated as 1X."),
    )
    build.add_argument("--output", required=True, type=Path, help="Output PON pickle path.")
    build.add_argument(
        "--stats-output",
        type=Path,
        help="Optional JSON or TSV file with PON build statistics.",
    )
    build.add_argument(
        "--pattern",
        default="*.tsv",
        help="Glob pattern for normal depth files inside --depth-dir.",
    )
    build.add_argument("--max-pca-components", default=3, type=int)

    call = subparsers.add_parser("call-cnv", help="Call one sample against a PON.")
    call.add_argument("--depth", required=True, type=Path, help="Sample depth TSV.")
    call.add_argument("--pon", required=True, type=Path, help="Serialized PON model.")
    call.add_argument(
        "--sample",
        help=(
            "Sample identifier used in reports, plot titles, and default output filenames. "
            "Defaults to the depth-file stem."
        ),
    )
    call.add_argument(
        "--output",
        type=Path,
        help="Output JSON or TSV report. Defaults to <sample>.echo.json.",
    )
    call.add_argument(
        "--cnv-output",
        type=Path,
        help="Optional dedicated TSV containing gene, exon, and breakpoint calls.",
    )
    call.add_argument(
        "--gene-output",
        type=Path,
        help="Optional compact TSV containing one gene-level copy-number row per gene.",
    )
    call.add_argument(
        "--plot",
        type=Path,
        help="Optional PNG/PDF/SVG path for a CNV, z-score, and PDS signal plot.",
    )

    return parser


def run_build_pon(args: argparse.Namespace) -> None:
    """Run the build-pon subcommand."""

    samples = _build_pon_samples(args)
    builder = PanelOfNormalsBuilder(
        args.bed,
        tiling_path=args.tiling,
        max_pca_components=args.max_pca_components,
    )
    pon = builder.build_from_samples(samples)
    pon.save(args.output)
    if args.stats_output is not None:
        pon.write_stats(args.stats_output)
    print(
        f"Wrote PON to {args.output} with {len(samples)} normals, "
        f"{len(pon.region_stats)} regions, modality={pon.modality}, "
        f"pds_mode={pon.pds_mode}, tiled_genes={len(pon.tiling_factors)}"
    )


def _build_pon_samples(args: argparse.Namespace) -> list[DepthSample]:
    if args.manifest is not None:
        return load_depth_manifest(args.manifest)
    if args.depth_dir is None:
        raise ValueError("build-pon requires either --manifest or --depth-dir")
    depth_paths = sorted(args.depth_dir.glob(args.pattern))
    if not depth_paths:
        raise FileNotFoundError(f"No depth files matching {args.pattern!r} in {args.depth_dir}")
    return [DepthSample(sample=path.stem, depth_file=path) for path in depth_paths]


def run_call_cnv(args: argparse.Namespace) -> None:
    """Run the call-cnv subcommand."""

    pon = PONModel.load(args.pon)
    caller = CNVCaller(pon)
    sample = args.sample or args.depth.stem
    report_path = args.output or Path(f"{sample}.echo.json")
    caller.write_outputs(
        args.depth,
        report_path=report_path,
        cnv_output_path=args.cnv_output,
        gene_output_path=args.gene_output,
        plot_path=args.plot,
        sample=sample,
    )
    print(f"Wrote ECHO report to {report_path}")
    if args.cnv_output is not None:
        print(f"Wrote CNV calls to {args.cnv_output}")
    if args.gene_output is not None:
        print(f"Wrote gene copy numbers to {args.gene_output}")
    if args.plot is not None:
        print(f"Wrote CNV plot to {args.plot}")


def main(argv: list[str] | None = None) -> int:
    """CLI entry point."""

    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "build-pon":
        run_build_pon(args)
    elif args.command == "call-cnv":
        run_call_cnv(args)
    else:
        parser.error(f"Unknown command: {args.command}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
