#!/usr/bin/env python3
"""
Generate and run DockThor TOML configurations for a qualification set.

The qualification set directory is expected to contain:

    config/<pdb>.grid_center.txt   three floats (x y z)
    config/<pdb>.grid_size.txt     three floats (x y z), in Angstroms
    config/<pdb>.rstep.txt         one float (grid step)
    ligand/<pdb>_ligand.top        ligand topology
    receptor/<pdb>.in              receptor topology
    cofactor/<pdb>/                optional, cofactor topologies

`--generate` writes one `<pdb>.toml` per system, in the root of the set or in
`--output-toml-path`. Input paths inside the TOMLs are relative to the set root.
`--dockthor-path <dockthor>` runs DockThor, from the set root, on every `<pdb>.toml`
found there, always with DockThor's `--success-rate` (added when `--dockthor-args` lacks it), then writes
qualification_summary.csv, the consolidation of every <ligand>_docked_success_rate.csv,
and success_rate.txt, the % of complexes whose topK_RMSD <= cutoff for every topK (and whose
overall_min_RMSD <= cutoff, column overall) and every cutoff 0.5, 1.0, ... up to `--rmsd-cutoff`. Both can be given
in the same call.
`--generate-summary <qualification_summary.csv>` only writes success_rate.txt, in the current directory,
from an existing summary (e.g. one joined by hand from interrupted runs); it docks nothing.
`--collect <dir>` writes both qualification_summary.csv and success_rate.txt from every
<ligand>_docked_success_rate.csv under <dir> (e.g. the reports of parallel CI jobs gathered in one folder);
it docks nothing either.
`--git-url <repo>` clones the set from a git repository (GitHub or any git URL)
and runs the same workflow on the clone. `--config-dir <dir>` runs TOMLs that already
exist (e.g. the `toml/` folder of a dataset repository) instead of generated ones.

Examples:
    python3 scripts/qualification_set.py ~/vcoreset-fix --generate -o ~/results/vcoreset
    python3 scripts/qualification_set.py ~/vcoreset-fix --dockthor-path build/bin/dockthor --mpi --np 8
    python3 scripts/qualification_set.py ~/vcoreset-fix --generate -o out --set ga-nrun=30 \\
        --dockthor-path build/bin/dockthor --systems 1bcu 10gs
    python3 scripts/qualification_set.py ~/vcoreset-fix --generate -o steadystate --output-toml-path tomls/ss \\
        --set success-rate=true --dockthor-path build/bin/dockthor
    python3 scripts/qualification_set.py --git-url https://github.com/<org>/<dataset>.git --generate -o steadystate \\
        --set success-rate=true --dockthor-path build/bin/dockthor --mpi --np 2
    python3 scripts/qualification_set.py --git-url git@github.com:<org>/<dataset>.git --config-dir toml \\
        -o steadystate --dockthor-path build/bin/dockthor \\
        --dockthor-args='--ga-nrun 24 --success-rate --success-rate-by-run 2.0'
    python3 scripts/qualification_set.py --config-dir ~/dockthor_coreset/toml -o steadystate \\
        --dockthor-path build/bin/dockthor --rmsd-cutoff 3
    python3 scripts/qualification_set.py --generate-summary steadystate/qualification_summary.csv --rmsd-cutoff 3
    python3 scripts/qualification_set.py --collect ci_collect/steadystate_symmetric --config-dir toml
"""

import argparse
import sys
from pathlib import Path

# The library lives next to this script (scripts/coreset/); scripts/ci/badge.py imports this module by path.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from coreset import docking, summary  # noqa: E402
from coreset.dataset import build_config, clone_dataset, discover_systems, infer_dataset  # noqa: E402
from coreset.summary import read_summary  # noqa: E402, F401  (scripts/ci/badge.py: qs.read_summary)
from coreset.tomlfile import parse_set_option, write_toml  # noqa: E402


def generate(args, systems: list[str]) -> int:
    output_dir = args.output_dir.expanduser()
    extra = dict(args.set)
    failures = 0
    written = 0
    args.toml_dir.mkdir(parents=True, exist_ok=True)
    for pdb in systems:
        toml_path = args.toml_dir / f"{pdb}.toml"
        if toml_path.exists() and not args.overwrite:
            print(f"[skip] {toml_path.name} exists (use --overwrite)")
            continue
        try:
            config = build_config(args.dataset, pdb, output_dir, not args.no_cofactors)
        except (OSError, ValueError) as error:
            print(f"[error] {pdb}: {error}", file=sys.stderr)
            failures += 1
            continue
        config.update(extra)
        write_toml(toml_path, pdb, config)
        written += 1
    print(f"Generated {written} TOML file(s) in {args.toml_dir} ({failures} error(s)).")
    return failures


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate and run DockThor TOML configurations for a qualification set.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "dataset",
        type=Path,
        nargs="?",
        help="Qualification set root (contains config/, ligand/, receptor/, ...). "
        "Omit it with --git-url, or with --config-dir to infer it from the TOMLs.",
    )
    parser.add_argument("--systems", nargs="+", metavar="PDB", help="Only these PDB IDs (default: all in config/).")

    src = parser.add_argument_group("dataset from git")
    src.add_argument(
        "--git-url",
        metavar="URL",
        help="Clone the qualification set from this git repository (e.g. "
        "https://github.com/<org>/<repo>.git) and run the workflow on the clone.",
    )
    src.add_argument(
        "--git-ref", metavar="REF", help="Branch, tag or commit to check out (default: the default branch)."
    )
    src.add_argument(
        "--git-subdir", metavar="PATH", help="Set root inside the repository, if it is not the repository root."
    )
    src.add_argument(
        "--clone-dir",
        type=Path,
        metavar="DIR",
        help="Where to clone (default: ./<repo name>). An existing clone there is reused as is.",
    )
    src.add_argument(
        "--git-update",
        action="store_true",
        help="Update an existing clone to --git-ref (or the remote default branch) before running.",
    )

    gen = parser.add_argument_group("generation")
    gen.add_argument("--generate", action="store_true", help="Write the <pdb>.toml files.")
    gen.add_argument(
        "--output-toml-path",
        type=Path,
        metavar="DIR",
        help="Directory of the TOML files, written by --generate and read by --dockthor-path "
        "(default: the dataset root). Relative to the current directory. The ligand, receptor "
        "and cofactor paths inside the TOMLs stay relative to the dataset root.",
    )
    gen.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        help="DockThor output directory; DockThor writes each system to <output-dir>/<pdb>. "
        "Written into the TOMLs with --generate; with --dockthor-path it is also passed to DockThor "
        "as --output-dir, overriding the TOMLs. A relative path is relative to the dataset root.",
    )
    gen.add_argument(
        "--set",
        action="append",
        default=[],
        type=parse_set_option,
        metavar="KEY=VALUE",
        help="Extra DockThor option written to every TOML, e.g. --set ga-nrun=30 "
        "--set seed=7 --set 'ga-evolution-type=\"GENERATIONAL_OCL\"'. Repeatable.",
    )
    gen.add_argument("--no-cofactors", action="store_true", help="Do not add cofactors-path.")
    gen.add_argument("--overwrite", action="store_true", help="Overwrite existing TOML files.")

    exe = parser.add_argument_group("execution")
    exe.add_argument(
        "--dockthor-path", type=Path, metavar="DOCKTHOR", help="Path to the dockthor binary; runs every TOML."
    )
    exe.add_argument(
        "--config-dir",
        type=Path,
        metavar="DIR",
        help="Run the TOMLs that already exist in DIR (every <pdb>.toml, or --systems) instead of "
        "the generated ones. Their input paths must be relative to the dataset root. With a "
        "dataset or --git-url, a relative DIR is relative to the dataset root (e.g. 'toml'); "
        "alone, DIR is relative to the current directory and the dataset root is inferred "
        "from the TOMLs (DIR or the first parent holding their ligand and receptor).",
    )
    exe.add_argument("--mpi", action="store_true", default=False, help="Run DockThor with MPI.")
    exe.add_argument("--np", type=int, default=4, help="Number of MPI processes (with --mpi).")
    exe.add_argument("--mpirun", default="mpirun", help="MPI launcher command.")
    exe.add_argument(
        "--mpi-args", default="", help="Extra launcher arguments; use the = form, e.g. --mpi-args=--oversubscribe."
    )
    exe.add_argument(
        "--dockthor-args",
        default="",
        help="Extra command-line arguments passed to DockThor; use the = form, e.g. --dockthor-args='--seed 7'. "
        "--success-rate is always added.",
    )
    exe.add_argument(
        "--summary-csv",
        type=Path,
        metavar="FILE",
        help="Consolidation of every <pdb>/<ligand>_docked_success_rate.csv (DockThor --success-rate): "
        "complex, execution_time(s), num_top, overall_min_RMSD_E, overall_min_RMSD (energy and RMSD "
        "of the lowest-RMSD leader of the all-runs clustering), top1_E, top1_RMSD, ..., "
        "topN_E, topN_RMSD, one row "
        "per complex (default: qualification_summary.csv in the common output directory, else "
        "next to the TOMLs). success_rate.txt is written next to it.",
    )
    exe.add_argument(
        "--rmsd-cutoff",
        type=float,
        default=2.0,
        metavar="RMSD",
        help="Largest RMSD cutoff of success_rate.txt: for each cutoff 0.5, 1.0, ... up to RMSD "
        "(rounded up to 0.5) and each topK, the %% of complexes with topK_RMSD <= cutoff, and in the "
        "overall column with overall_min_RMSD <= cutoff.",
    )
    exe.add_argument(
        "--generate-summary",
        type=Path,
        metavar="CSV",
        help="Only write success_rate.txt, in the current directory, from the existing "
        "qualification_summary.csv CSV (with --rmsd-cutoff); nothing is generated or docked.",
    )
    exe.add_argument(
        "--collect",
        type=Path,
        metavar="DIR",
        help="Only consolidate: write qualification_summary.csv (or --summary-csv) and success_rate.txt "
        "from every <ligand>_docked_success_rate.csv under DIR, flat or in <pdb>/ folders (e.g. "
        "the reports of parallel runs gathered in one folder). With --config-dir, the complexes "
        "are counted against its TOMLs and the missing ones listed. Nothing is generated or docked.",
    )
    exe.add_argument(
        "--skip-done", action="store_true", help="Skip systems that already have a *_docked.csv in <output-dir>/<pdb>."
    )
    exe.add_argument("--stop-on-error", action="store_true", help="Stop at the first failing system.")
    exe.add_argument("--dry-run", action="store_true", help="Print the commands without running them.")
    exe.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Also print DockThor output to the terminal; it is always saved to " "<output-dir>/<pdb>/<pdb>_run.log.",
    )

    args = parser.parse_args(argv)
    if args.collect is not None:
        if (
            args.dataset is not None
            or args.git_url
            or args.generate
            or args.dockthor_path
            or args.generate_summary is not None
        ):
            parser.error(
                "--collect runs alone (with --config-dir, --summary-csv and --rmsd-cutoff): drop the "
                "dataset, --git-url, --generate, --dockthor-path and --generate-summary"
            )
        if args.rmsd_cutoff <= 0:
            parser.error("--rmsd-cutoff must be positive")
        args.parser = parser
        return args
    if args.generate_summary is not None:
        if args.dataset is not None or args.git_url or args.config_dir or args.generate or args.dockthor_path:
            parser.error(
                "--generate-summary runs alone: drop the dataset, --git-url, --config-dir, "
                "--generate and --dockthor-path"
            )
        if args.rmsd_cutoff <= 0:
            parser.error("--rmsd-cutoff must be positive")
        args.parser = parser
        return args
    if args.dataset is not None and args.git_url is not None:
        parser.error("give either the dataset directory or --git-url, not both")
    if args.dataset is None and args.git_url is None and args.config_dir is None:
        parser.error(
            "give the dataset directory, --git-url, or --config-dir with the TOMLs " "(or --generate-summary CSV)"
        )
    if args.dataset is None and args.git_url is None and args.generate:
        parser.error("--generate needs the dataset directory or --git-url")
    if args.git_url is None and (args.git_ref or args.git_subdir or args.clone_dir or args.git_update):
        parser.error("--git-ref, --git-subdir, --clone-dir and --git-update need --git-url")
    if not args.generate and args.dockthor_path is None and args.git_url is None:
        parser.error("nothing to do: pass --generate and/or --dockthor-path")
    if args.generate and args.output_dir is None:
        parser.error("--generate requires --output-dir")
    if args.config_dir is not None and args.dockthor_path is None:
        parser.error("--config-dir requires --dockthor-path")
    if args.rmsd_cutoff <= 0:
        parser.error("--rmsd-cutoff must be positive")
    args.parser = parser
    return args


def resolve_dataset(args) -> None:
    """Set args.dataset (cloning it first with --git-url, or inferring it from --config-dir) and args.toml_dir."""
    if args.dataset is None and args.git_url is None:
        toml_dir = args.config_dir.expanduser().resolve()
        if not toml_dir.is_dir():
            args.parser.error(f"--config-dir {toml_dir} is not a directory")
        dataset = infer_dataset(toml_dir)
        if dataset is None:
            args.parser.error(
                f"cannot find the dataset root of the TOMLs in {toml_dir}: no folder from there up "
                "holds their ligand and receptor; pass the dataset directory"
            )
        print(f"Dataset root (from the TOMLs): {dataset}")
        args.config_dir = toml_dir  # absolute, so it does not depend on the dataset root below
    else:
        dataset = clone_dataset(args) if args.git_url else args.dataset
    args.dataset = dataset.expanduser().resolve()
    args.toml_dir = args.output_toml_path.expanduser().resolve() if args.output_toml_path else args.dataset
    args.run_toml_dir = None
    if args.config_dir is not None:
        args.run_toml_dir = (args.dataset / args.config_dir.expanduser()).resolve()  # an absolute DIR stays as is
        if not args.run_toml_dir.is_dir():
            args.parser.error(f"--config-dir {args.run_toml_dir} is not a directory")
    needs_config = args.generate or (args.dockthor_path is not None and args.config_dir is None and not args.systems)
    if needs_config and not (args.dataset / "config").is_dir():
        args.parser.error(f"{args.dataset} has no config/ directory")


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.collect is not None:
        return summary.collect(args.collect, args.config_dir, args.summary_csv, args.rmsd_cutoff)
    if args.generate_summary is not None:
        return summary.generate_summary(args.generate_summary, args.rmsd_cutoff)
    resolve_dataset(args)
    errors = 0
    if args.generate:
        errors += generate(args, args.systems or discover_systems(args.dataset))
    if args.dockthor_path is not None:
        errors += docking.run(args, docking.run_systems(args))
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
