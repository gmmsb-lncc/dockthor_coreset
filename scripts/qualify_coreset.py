#!/usr/bin/env python3
"""Build a qualification set from the CASF-2013 core set with the DockThor Suite, then dock it.

The whole chain, from the files CASF-2013 distributes to the success rates:

    prepare   MMFFLigand on <pdb>_ligand.mol2 and on each metal ion, PDBThorBox on <pdb>_protein.pdb;
              everything they write lands in --output-dir, the set's files in --dataset
    config    grid_toml_generator.py            -> <dataset>/config/<pdb>.toml
    dock      DockThor, through qualification_set.py
                                                -> <dataset>/<--dock-output>/, qualification_summary.csv,
                                                   success_rate.txt

Every stage runs by default; --skip-mmffligand, --skip-pdbthorbox, --skip-config and --skip-dockthor leave one out
(skipping both MMFFLigand and PDBThorBox skips the prepare stage). What a skipped stage would have written is taken
from the set as it is.

Two directories, because the facilities write much more than the set carries:

    --output-dir  one folder per complex with **every** generated file -- the ligand topology and
                  <pdb>_ligand_mmff.mol2, the receptor topology with its charged PDB and MOL2, the chain
                  files, the residue and summary reports, resumo.out, the ion topologies in cofactors/,
                  and the console output of every call in logs/
    --dataset     the set itself (README 1.4): <pdb>/ with receptors/, ligands/ and cofactors/, then config/

A set directory holds one qualification set, so several can live side by side in the repository
(`v2013-coreset/`, ...). Each stage can be skipped, so a failed one is repeated without redoing the others.

    python3 scripts/qualify_coreset.py ~/Downloads/CASF-2013/coreset \
        --output-dir ~/prepared_v2013 --dataset v2013-coreset --suite-bin $DOCKTHOR/build/bin --skip-dockthor
    python3 scripts/qualify_coreset.py ~/Downloads/CASF-2013/coreset --dataset v2013-coreset \
        --skip-mmffligand --skip-pdbthorbox --skip-config --dockthor-path $DOCKTHOR/build/bin/dockthor \
        --dockthor-args='--ga-nrun 24 --success-rate --success-rate-by-run 3 --symmetric-RMSD'

The ligands are the deposited poses (`<pdb>_ligand.mol2`), never `<pdb>_ligand_opt.mol2`: the RMSD DockThor
reports is measured against the input topology, so the reference has to be the crystal pose (README 1.2).
"""

import argparse
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from coreset.dataset import ligand_path, receptor_path  # noqa: E402
from coreset.prepare import (  # noqa: E402
    Prepared,
    StepFailed,
    Tools,
    discover_casf,
    prepare_complex,
    write_index,
)

SCRIPTS = Path(__file__).resolve().parent


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a qualification set from the CASF-2013 core set and dock it.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("casf", type=Path, help="CASF-2013 core-set directory (one folder per complex)")
    parser.add_argument(
        "-d",
        "--dataset",
        type=Path,
        required=True,
        help="Directory of the set to build: <pdb>/receptors/, <pdb>/ligands/, <pdb>/cofactors/, config/",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Directory for every file the facilities write, one folder per complex "
        "(default: <dataset>_prepared beside the set). Only the prepare stage uses it",
    )
    parser.add_argument("--complex", nargs="+", metavar="PDB", help="Only these complexes (default: all)")
    skip = parser.add_argument_group("stages (all run by default)")
    skip.add_argument(
        "--skip-mmffligand",
        action="store_true",
        help="Do not run MMFFLigand: keep the set's ligands and cofactors (prepare stage)",
    )
    skip.add_argument(
        "--skip-pdbthorbox", action="store_true", help="Do not run PDBThorBox: keep the set's receptors (prepare stage)"
    )
    skip.add_argument(
        "--skip-config", action="store_true", help="Do not run grid_toml_generator.py: keep the set's config/"
    )
    skip.add_argument("--skip-dockthor", action="store_true", help="Do not dock (dock stage)")
    parser.add_argument(
        "--suite-bin",
        type=Path,
        default=Path.home() / "dockthor-suite" / "build" / "bin",
        help="Directory holding mmffligand, pdbthorbox and dockthor",
    )
    parser.add_argument(
        "--extra-cofactors",
        type=Path,
        help="Another set whose <pdb>/cofactors/ is copied for a complex whose CASF protein file holds no "
        "ion (the Hg of 1n1m, the Cd of 3ao4 and the Ni of 3cj2, README 1.6)",
    )
    parser.add_argument("--jobs", type=int, default=4, help="Complexes prepared at the same time")
    parser.add_argument("--padding", type=float, default=8.0, help="Å added to each box edge (config stage)")
    parser.add_argument("--box-shape", choices=("cube", "ligand"), default="cube", help="Box shape (config stage)")
    parser.add_argument("--rstep", type=float, help="Grid step written into the TOMLs (config stage)")
    parser.add_argument(
        "--dockthor-path", type=Path, help="The dockthor binary (dock stage; default: <--suite-bin>/dockthor)"
    )
    parser.add_argument("--dock-output", default="steadystate", help="-o of the dock stage, under the set")
    parser.add_argument("--dockthor-args", default="", help="Extra DockThor options (dock stage)")
    parser.add_argument("--mpi", action="store_true", help="Launch DockThor with mpirun (dock stage)")
    parser.add_argument("--np", type=int, default=4, help="MPI processes per complex (dock stage)")
    parser.add_argument("--overwrite", action="store_true", help="Rebuild files that already exist")
    args = parser.parse_args(argv)

    args.suite_bin = args.suite_bin.expanduser().resolve()
    args.dockthor_path = (
        args.dockthor_path.expanduser().resolve() if args.dockthor_path is not None else args.suite_bin / "dockthor"
    )
    args.stages = []
    if not (args.skip_mmffligand and args.skip_pdbthorbox):
        args.stages.append("prepare")
    if not args.skip_config:
        args.stages.append("config")
    if not args.skip_dockthor:
        args.stages.append("dock")
    if not args.stages:
        parser.error("every stage is skipped, nothing to do")

    # The binaries are checked before anything runs, so a missing one does not stop the chain hours later.
    needed = []
    if not args.skip_mmffligand:
        needed.append((args.suite_bin / "mmffligand", "--suite-bin", "--skip-mmffligand"))
    if not args.skip_pdbthorbox:
        needed.append((args.suite_bin / "pdbthorbox", "--suite-bin", "--skip-pdbthorbox"))
    if not args.skip_dockthor:
        needed.append((args.dockthor_path, "--dockthor-path", "--skip-dockthor"))
    for binary, option, skip_option in needed:
        if not binary.is_file():
            parser.error(f"no {binary}: set {option}, or give {skip_option}")
    return args


def already_prepared(args: argparse.Namespace, pdb: str) -> bool:
    """The set holds what the facilities that run would write: the ligand (MMFFLigand), the receptor (PDBThorBox)."""
    ligand = args.skip_mmffligand or (args.dataset / ligand_path(pdb)).is_file()
    receptor = args.skip_pdbthorbox or (args.dataset / receptor_path(pdb)).is_file()
    return ligand and receptor


def step_prepare(args: argparse.Namespace, systems: list[str]) -> int:
    """MMFFLigand (ligand and cofactors) and PDBThorBox (receptor) of every complex, unless skipped. The number that
    failed."""
    tools = Tools.find(args.suite_bin, mmffligand=not args.skip_mmffligand, pdbthorbox=not args.skip_pdbthorbox)
    facilities = " and ".join(
        name for name, binary in (("MMFFLigand", tools.mmffligand), ("PDBThorBox", tools.pdbthorbox)) if binary
    )
    print(f"\n== prepare: {len(systems)} complex(es), {facilities} of {args.suite_bin}")
    print(f"   generated files -> {args.output_dir}")
    print(f"   set             -> {args.dataset}")

    todo = systems
    if not args.overwrite:
        todo = [pdb for pdb in systems if not already_prepared(args, pdb)]
        if len(todo) < len(systems):
            print(f"   {len(systems) - len(todo)} already prepared, kept (use --overwrite to rebuild)")

    done: list[Prepared] = []
    failed: list[tuple[str, str]] = []

    def one(pdb: str) -> tuple[str, Prepared | str]:
        try:
            return pdb, prepare_complex(tools, args.casf, pdb, args.output_dir, args.dataset, args.extra_cofactors)
        except (StepFailed, OSError, ValueError) as error:
            return pdb, str(error)

    started = time.monotonic()
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        for number, (pdb, result) in enumerate(pool.map(one, todo), start=1):
            if isinstance(result, Prepared):
                done.append(result)
                parts = []
                if result.ligand_atoms is not None:
                    extra = f" (+{result.carried_cofactors} carried)" if result.carried_cofactors else ""
                    parts += [f"ligand {result.ligand_atoms} atoms", f"{result.cofactors} cofactor(s){extra}"]
                if result.receptor_atoms is not None:
                    parts.append(f"receptor {result.receptor_atoms}")
                    if result.waters_dropped:
                        parts.append(f"{result.waters_dropped} water(s) dropped")
                print(f"[{number}/{len(todo)}] {pdb}: " + ", ".join(parts))
            else:
                failed.append((pdb, result))
                print(f"[{number}/{len(todo)}] {pdb}: [error] {result}", file=sys.stderr)

    print(f"   prepared {len(done)}, failed {len(failed)}, in {time.monotonic() - started:.1f}s")
    if written := write_index(args.casf, args.dataset):
        print(f"   {written} affinities written to {args.dataset / 'index' / 'affinity.csv'}")
    if done:
        if not args.skip_mmffligand:
            with_cofactors = [p for p in done if p.cofactors or p.carried_cofactors]
            total = sum(p.cofactors + p.carried_cofactors for p in done)
            print(f"   {len(with_cofactors)} complex(es) with cofactors, {total} topologies")
        print(f"   {sum(p.generated_files for p in done)} generated file(s) in {args.output_dir}")
    if failed:
        print("   failed: " + " ".join(pdb for pdb, _ in failed), file=sys.stderr)
    return len(failed)


def step_config(args: argparse.Namespace, systems: list[str]) -> int:
    """grid_toml_generator.py over the prepared ligands -> config/<pdb>.toml."""
    command = [
        sys.executable,
        str(SCRIPTS / "grid_toml_generator.py"),
        str(args.dataset),
        "--padding",
        str(args.padding),
        "--box-shape",
        args.box_shape,
    ]
    if args.rstep is not None:
        command += ["--rstep", str(args.rstep)]
    if args.complex:
        command += ["--complex", *systems]
    if args.overwrite:
        command.append("--overwrite")
    print(f"\n== config: {' '.join(command[1:])}")
    return subprocess.run(command).returncode


def step_dock(args: argparse.Namespace, systems: list[str]) -> int:
    """qualification_set.py on the generated TOMLs."""
    command = [
        sys.executable,
        str(SCRIPTS / "qualification_set.py"),
        str(args.dataset),
        "--config-dir",
        "config",
        "-o",
        args.dock_output,
        "--dockthor-path",
        str(args.dockthor_path),
        "--skip-done",
    ]
    if args.mpi:
        command += ["--mpi", "--np", str(args.np)]
    if args.complex:
        command += ["--complex", *systems]
    if args.dockthor_args:
        command.append(f"--dockthor-args={args.dockthor_args}")
    print(f"\n== dock: {' '.join(command[1:])}")
    return subprocess.run(command).returncode


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.casf = args.casf.expanduser().resolve()
    args.dataset = args.dataset.expanduser().resolve()
    # The generated files default to a folder beside the set, named after it.
    args.output_dir = (
        args.output_dir.expanduser().resolve()
        if args.output_dir is not None
        else args.dataset.with_name(args.dataset.name + "_prepared")
    )
    if args.extra_cofactors is not None:
        args.extra_cofactors = args.extra_cofactors.expanduser().resolve()
    if not args.casf.is_dir():
        sys.exit(f"[error] {args.casf} is not a directory")

    available = discover_casf(args.casf)
    if not available:
        sys.exit(f"[error] {args.casf} holds no <pdb>/<pdb>_protein.pdb + <pdb>_ligand.mol2")
    systems = args.complex or available
    unknown = sorted(set(systems) - set(available))
    if unknown:
        sys.exit(f"[error] not in {args.casf}: {' '.join(unknown)}")

    print(f"CASF core set: {args.casf} ({len(available)} complexes, {len(systems)} selected)")
    print(f"Set to build:  {args.dataset}")
    args.dataset.mkdir(parents=True, exist_ok=True)
    if "prepare" in args.stages:
        args.output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Stages:        {', '.join(args.stages)}")

    status = 0
    for stage in args.stages:
        failures = {"prepare": step_prepare, "config": step_config, "dock": step_dock}[stage](args, systems)
        if failures:
            status = 1
            if stage != "dock":
                print(f"[error] the {stage} stage failed; later stages are skipped", file=sys.stderr)
                break
    return status


if __name__ == "__main__":
    sys.exit(main())
