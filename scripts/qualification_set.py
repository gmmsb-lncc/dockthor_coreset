#!/usr/bin/env python3
"""
Dock a qualification set with DockThor and report its success rates.

The set directory holds, per complex (README 1.4):

    <pdb>/receptors/<pdb>.in       receptor topology
    <pdb>/receptors/prep/*_prep.pdb charged PDB, the receptor DockTScore reads
    <pdb>/ligands/<pdb>_ligand.top ligand topology
    <pdb>/cofactors/*.top          optional, cofactor topologies (and *.pdb, read by DockTScore)
    index/affinity.csv             the measured affinities
    config/<pdb>.toml              the inputs and the box, written by scripts/grid_toml_generator.py

DockThor runs from the set root on every <pdb>.toml of --config-dir (default config/, or --complex), always with
--success-rate (added when --dockthor-args lacks it). After every complex the script rewrites
qualification_summary.csv, the consolidation of every <ligand>_docked_success_rate.csv, and success_rate.txt, the %
of complexes whose topK_RMSD <= cutoff for every topK (and whose overall_min_RMSD <= cutoff, column overall) and every
cutoff 0.5, 1.0, ... up to --rmsd-cutoff. Reports docked elsewhere are merged with scripts/summarize_results.py.

The summary always carries the measured affinity of each complex (log_affinity, delta_g_kcal_mol, from
index/affinity.csv). With --docktscore-path, DockTScore also rescores the best-scored pose of every complex right
after it docks, once per model of --docktscore-models: the scores go in the summary as docktscore_<model> columns,
and scoring_power.txt and scoring_power.csv beside success_rate.txt report the Pearson r, R^2, Spearman rho and the
two RMSE of every model against the measured affinity. A finished run is rescored afterwards, without docking it
again, with scripts/rescore_results.py.

The set is a directory, a git repository (--git-url, cloned or reused), or inferred from --config-dir alone (the
first folder from there up that holds the TOMLs' ligand and receptor).

Examples:
    python3 scripts/qualification_set.py v2013-coreset -o steadystate --dockthor-path build/bin/dockthor --mpi --np 8
    python3 scripts/qualification_set.py v2013-coreset --complex 1bcu 10gs -o check \\
        --dockthor-path build/bin/dockthor --dockthor-args='--ga-nrun 4 --seed 7'
    python3 scripts/qualification_set.py --config-dir v2013-coreset/config -o steadystate \\
        --dockthor-path build/bin/dockthor --rmsd-cutoff 3
    python3 scripts/qualification_set.py v2013-coreset -o scored --dockthor-path build/bin/dockthor \\
        --docktscore-path build/bin/docktscore            # dock, rescore, and report the scoring power
    python3 scripts/qualification_set.py --git-url https://github.com/<org>/<dataset>.git --git-subdir v2013-coreset \\
        -o steadystate --dockthor-path build/bin/dockthor \\
        --dockthor-args='--ga-nrun 24 --success-rate --success-rate-by-run 2.0'

The TOMLs are generated first, with scripts/grid_toml_generator.py <set>.
"""

import argparse
import shlex
import sys
from pathlib import Path

# The library lives next to this script (scripts/coreset/).
sys.path.insert(0, str(Path(__file__).resolve().parent))

from coreset import docking, scoring  # noqa: E402
from coreset.dataset import clone_dataset, infer_dataset  # noqa: E402

DEFAULT_CONFIG_DIR = Path("config")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Dock a qualification set with DockThor and report its success rates.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "dataset",
        type=Path,
        nargs="?",
        help="Qualification set root (one <pdb>/ directory per complex, and config/). "
        "Omit it with --git-url, or with --config-dir to infer it from the TOMLs.",
    )
    parser.add_argument(
        "--complex", nargs="+", metavar="PDB", help="Only these PDB IDs (default: every TOML of --config-dir)."
    )

    src = parser.add_argument_group("dataset from git")
    src.add_argument(
        "--git-url",
        metavar="URL",
        help="Clone the qualification set from this git repository (e.g. https://github.com/<org>/<repo>.git) "
        "and dock the clone.",
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

    exe = parser.add_argument_group("docking")
    exe.add_argument(
        "--dockthor-path", type=Path, required=True, metavar="DOCKTHOR", help="Path to the dockthor binary."
    )
    exe.add_argument(
        "--config-dir",
        type=Path,
        metavar="DIR",
        help=f"The TOMLs to run (every <pdb>.toml, or --complex). With a dataset or --git-url, a relative DIR is "
        f"relative to the dataset root (default: {DEFAULT_CONFIG_DIR}); alone, DIR is relative to the current "
        "directory and the dataset root is inferred from the TOMLs.",
    )
    exe.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        help="DockThor output directory, passed as --output-dir (overriding the TOMLs' output-dir); DockThor writes "
        "each complex to <output-dir>/<pdb>. A relative path is relative to the dataset root.",
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
        help="Consolidation of every <pdb>/<ligand>_docked_success_rate.csv (DockThor --success-rate), one row per "
        "complex (default: qualification_summary.csv in the common output directory, else next to the TOMLs). "
        "success_rate.txt is written next to it.",
    )
    exe.add_argument(
        "--rmsd-cutoff",
        type=float,
        default=3.0,
        metavar="RMSD",
        help="Largest RMSD cutoff of success_rate.txt: for each cutoff 0.5, 1.0, ... up to RMSD (rounded up to 0.5) "
        "and each topK, the %% of complexes with topK_RMSD <= cutoff, and in the overall column with "
        "overall_min_RMSD <= cutoff. Default 3.0.",
    )
    score = parser.add_argument_group("rescoring with DockTScore")
    score.add_argument(
        "--docktscore-path",
        type=Path,
        metavar="DOCKTSCORE",
        help="Path to the docktscore binary. Given, the best-scored pose of every complex is rescored right after "
        "it docks (<output-dir>/<pdb>/<pdb>_top1.mol2 and its .score), the score goes in the summary as the "
        "docktscore_<model> columns, and scoring_power.txt and scoring_power.csv are written beside "
        "success_rate.txt.",
    )
    score.add_argument(
        "--docktscore-models",
        "--docktscore-model",
        nargs="+",
        default=scoring.DEFAULT_MODELS,
        metavar="MODEL",
        help="DockTScore coefficients to score with, one summary column each (docktscore_<model>): general, "
        "general2764, general_dock, kinase, ppi, protease. Every default model scores a pose in about 0.1 s.",
    )
    score.add_argument(
        "--docktscore-args",
        default="",
        help="Extra command-line arguments passed to DockTScore; use the = form, e.g. --docktscore-args='-v'.",
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
        help="Also print DockThor output to the terminal; it is always saved to <output-dir>/<pdb>/<pdb>_run.log.",
    )

    args = parser.parse_args(argv)
    if args.dataset is not None and args.git_url is not None:
        parser.error("give either the dataset directory or --git-url, not both")
    if args.dataset is None and args.git_url is None and args.config_dir is None:
        parser.error("give the dataset directory, --git-url, or --config-dir with the TOMLs")
    if args.git_url is None and (args.git_ref or args.git_subdir or args.clone_dir or args.git_update):
        parser.error("--git-ref, --git-subdir, --clone-dir and --git-update need --git-url")
    if args.rmsd_cutoff <= 0:
        parser.error("--rmsd-cutoff must be positive")
    if args.docktscore_path is None and args.docktscore_args:
        parser.error("--docktscore-args needs --docktscore-path")
    args.docktscore_args = shlex.split(args.docktscore_args)
    if args.docktscore_path is not None:
        args.docktscore_path = args.docktscore_path.expanduser().resolve()
        if not args.docktscore_path.is_file():
            parser.error(f"--docktscore-path {args.docktscore_path} is not a file")
    args.parser = parser
    return args


def resolve_dataset(args) -> None:
    """Set args.dataset (cloning it first with --git-url, or inferring it from --config-dir) and args.run_toml_dir,
    the folder of the TOMLs to run."""
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
    config_dir = args.config_dir.expanduser() if args.config_dir is not None else DEFAULT_CONFIG_DIR
    args.run_toml_dir = (args.dataset / config_dir).resolve()  # an absolute DIR stays as is
    if not args.run_toml_dir.is_dir():
        args.parser.error(
            f"no TOML folder {args.run_toml_dir}; generate it with scripts/grid_toml_generator.py {args.dataset}"
        )


def main(argv=None) -> int:
    args = parse_args(argv)
    resolve_dataset(args)
    return 1 if docking.run(args, docking.run_systems(args)) else 0


if __name__ == "__main__":
    sys.exit(main())
