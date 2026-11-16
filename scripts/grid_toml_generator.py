#!/usr/bin/env python3
"""
Write a DockThor TOML for every complex of a set, with a docking box fitted to its ligand.

A complex is a directory of the set root holding its ligand (README 1.4):

    <complex>/receptors/<complex>.in            receptor topology
    <complex>/ligands/<complex>_ligand.top      ligand topology (or MMFFLigand's <complex>_ligand_mmff.cif /
                                                .mmcif, gzipped too: the formats DockThor reads as a ligand)
    <complex>/cofactors/*.top                   optional cofactor topologies

For each ligand, the box is:

    grid-center  the middle of the ligand's bounding box (all atoms, hydrogens included)
    grid-size    the bounding box's extent on each axis + --padding (half on each side);
                 --box-shape cube (default): every edge the longest one; --box-shape ligand: each axis its own
    grid-rstep   --rstep, written only when given (else DockThor's default, 0.25 Å, also used to count the points)

The values are float32, as DockThor reads them, written as the shortest decimal that reads back the same. The grid
points are counted as DockThor does (int(size / rstep) + 1 per axis, in float32); a grid over DockThor's default
limit of 1,000,000 points gets `grid-max-points` with that count, so DockThor runs it.

`<complex>.toml` is written to --output-dir with output-dir (--run-dir), ligand, receptor, cofactors-path (when
<complex>/cofactors/ holds .top files) and the grid. The ligand, receptor and cofactors-path are absolute, so DockThor
finds them from any working directory; generate the TOMLs on the machine (or in the container) that docks. The
output-dir stays relative, to the directory DockThor runs in (the set root with scripts/qualification_set.py, whose
-o overrides it). An existing TOML is kept unless --overwrite.

A file named <complex>_ligand* that is not a ligand topology is reported as `Skipped unrecognized <FILE>`;
MMFFLigand's deprecated copies (new_*, *_new.*) are ignored. When a ligand is there in several formats, the first of
.top, .cif, .mmcif, .cif.gz, .mmcif.gz is used.

Examples:
    python3 scripts/grid_toml_generator.py v2013-coreset                 # every complex, TOMLs in config/
    python3 scripts/grid_toml_generator.py v2013-coreset --padding 10 --box-shape ligand --complex 1bcu 10gs
    python3 scripts/grid_toml_generator.py v2013-coreset --output-dir toml_new --rstep 0.3 --overwrite
"""

import argparse
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from coreset.dataset import LIGAND_TAG, LIGANDS_DIR, cofactors_dir, complex_dirs, receptor_path  # noqa: E402
from coreset.gridbox import (  # noqa: E402
    DEFAULT_PADDING,
    DEFAULT_RSTEP,
    DOCKTHOR_MAX_GRID_POINTS,
    BoxShape,
    GridBox,
    grid_box,
)
from coreset.molecules import Ligand, find_ligands  # noqa: E402
from coreset.tomlfile import write_toml  # noqa: E402


def absolute(path: Path) -> str:
    """A topology path as written in a TOML: absolute, symbolic links resolved."""
    return path.resolve().as_posix()


def relative_to_root(path: Path, root: Path) -> str:
    """A path as printed in the messages: relative to the set root when it is inside it, else absolute."""
    try:
        return path.resolve().relative_to(root).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def toml_config(
    box: GridBox, ligand: str, receptor: str, cofactors: str | None, run_dir: str, write_rstep: bool
) -> tuple[dict, dict]:
    """The keys of the TOML, in the order of toml/, and the comment above grid-max-points when it is needed.

    write_rstep: write grid-rstep (when --rstep was given); without it DockThor uses its default step."""
    config = {"output-dir": run_dir, "ligand": ligand, "receptor": receptor}
    if cofactors is not None:
        config["cofactors-path"] = cofactors
    config["grid-size"] = list(box.size)
    config["grid-center"] = list(box.center)
    if write_rstep:
        config["grid-rstep"] = box.rstep
    comments = {}
    if box.total_points > DOCKTHOR_MAX_GRID_POINTS:
        x, y, z = box.points
        limit = DOCKTHOR_MAX_GRID_POINTS
        config["grid-max-points"] = box.total_points
        comments["grid-max-points"] = (
            f"grid points: {box.total_points} ({x} x {y} x {z}), over DockThor's default limit {limit}"
        )
    return config, comments


def vector_text(values, digits: int = 3) -> str:
    return " x ".join(f"{v:.{digits}f}" for v in values)


def positive(text: str) -> float:
    value = float(text)
    if not value > 0:
        raise argparse.ArgumentTypeError(f"must be positive, got {text}")
    return value


def non_negative(text: str) -> float:
    value = float(text)
    if not value >= 0:
        raise argparse.ArgumentTypeError(f"must be zero or positive, got {text}")
    return value


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Write a DockThor TOML for every ligand, with a docking box fitted to the ligand.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "dataset",
        type=Path,
        nargs="?",
        default=Path("."),
        help="Set root, one <complex>/ directory per complex.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("config"),
        help="Where to write <complex>.toml, relative to the set root.",
    )
    parser.add_argument(
        "--run-dir", default="vcoreset_run", help="DockThor's output directory, the output-dir key of the TOMLs."
    )
    parser.add_argument(
        "--padding",
        type=non_negative,
        default=DEFAULT_PADDING,
        help="Å added to each edge of the ligand's bounding box (half on each side).",
    )
    parser.add_argument(
        "--box-shape",
        type=BoxShape,
        choices=list(BoxShape),
        default=BoxShape.CUBE,
        help="cube: every edge the longest one; ligand: each edge the ligand's extent on that axis + padding.",
    )
    parser.add_argument(
        "--rstep",
        type=positive,
        help=f"Grid step (Å), written as grid-rstep only when given; without it DockThor's default {DEFAULT_RSTEP} "
        "is used, for the point count too.",
    )
    parser.add_argument("--complex", nargs="+", metavar="COMPLEX", help="Only these complexes (default: every ligand).")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing TOML files.")
    return parser.parse_args(argv)


def complex_ligand(root: Path, name: str) -> tuple[Ligand | None, list[Path]]:
    """The ligand of complex `name`, <name>/ligands/<name>_ligand in the first format DockThor reads, and the files
    named like it that are not ligand topologies."""
    ligands, skipped = find_ligands(root / name / LIGANDS_DIR)
    wanted = name + LIGAND_TAG
    ligand = next((ligand for ligand in ligands if ligand.name == wanted), None)
    return ligand, [path for path in skipped if path.name.startswith(wanted)]


def main(argv=None) -> int:
    args = parse_args(argv)
    root = args.dataset.expanduser().resolve()
    output_dir = root / args.output_dir.expanduser()
    if not root.is_dir():
        print(f"[error] {root} is not a directory", file=sys.stderr)
        return 1

    names = complex_dirs(root)
    errors = 0
    if args.complex:
        if unknown := [name for name in args.complex if name not in names]:
            print(f"[error] no <complex>/ligands/<complex>_ligand of {' '.join(unknown)} in {root}", file=sys.stderr)
            errors += len(unknown)
        names = [name for name in args.complex if name in names]

    by_complex = {}
    for name in names:
        ligand, skipped = complex_ligand(root, name)
        for path in skipped:
            print(f"Skipped unrecognized {relative_to_root(path, root)}", file=sys.stderr)
        if ligand is None:
            print(f"[error] {name}: no ligand topology DockThor reads in {name}/{LIGANDS_DIR}/", file=sys.stderr)
            errors += 1
            continue
        by_complex[name] = ligand

    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"{'Complex':<10} {'Bounding box (Å)':<26} {'Grid size (Å)':<26} {'Grid points':>12}")
    written = 0
    totals = []
    over_limit = []
    for name, ligand in by_complex.items():
        box = grid_box(ligand.coordinates, args.padding, args.box_shape, args.rstep or DEFAULT_RSTEP)
        totals.append(box.total_points)
        warning = ""
        if box.total_points > DOCKTHOR_MAX_GRID_POINTS:
            over_limit.append(name)
            warning = f"  WARNING: over {DOCKTHOR_MAX_GRID_POINTS} points, grid-max-points written"
        print(f"{name:<10} {vector_text(box.extent):<26} {vector_text(box.size):<26} {box.total_points:>12}{warning}")

        toml_path = output_dir / f"{name}.toml"
        if toml_path.exists() and not args.overwrite:
            print(f"[skip] {relative_to_root(toml_path, root)} exists (use --overwrite)")
            continue
        receptor = root / receptor_path(name)
        if not receptor.is_file():
            print(f"[error] {name}: missing {receptor}", file=sys.stderr)
            errors += 1
            continue
        cofactors = cofactors_dir(root, name)
        config, comments = toml_config(
            box,
            absolute(ligand.path),
            absolute(receptor),
            absolute(root / cofactors) if cofactors is not None else None,
            args.run_dir,
            write_rstep=args.rstep is not None,
        )
        write_toml(toml_path, name, config, comments)
        written += 1

    points = (
        f"; grid points min {min(totals)}, median {statistics.median(totals):.0f}, max {max(totals)}" if totals else ""
    )
    print(
        f"Wrote {written} TOML file(s) to {relative_to_root(output_dir, root)} for {len(by_complex)} ligand(s), "
        f"{len(over_limit)} over {DOCKTHOR_MAX_GRID_POINTS} points, {errors} error(s){points}."
    )
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
