#!/usr/bin/env python3
"""
Generate and run DockThor TOML configurations for a qualification set.

Part of the DockThor coreset repository (https://github.com/gmmsb-lncc/dockthor_coreset);
see its README.md for the full reference. Run from a clone, the dataset defaults to that clone.

The qualification set directory is expected to contain:

    config/<pdb>.grid_center.txt   three floats (x y z)
    config/<pdb>.grid_size.txt     three floats (x y z), in Angstroms
    config/<pdb>.rstep.txt         one float (grid step)
    ligand/<pdb>_ligand.top        ligand topology
    receptor/<pdb>.in              receptor topology
    cofactor/<pdb>/                optional, cofactor topologies

`--generate` writes one `<pdb>.toml` per system, in the root of the set or in
`--output-toml-path`. Input paths inside the TOMLs are relative to the set root.
`--run-path <dockthor>` runs DockThor, from the set root, on every `<pdb>.toml`
found there, then writes a summary CSV with the execution time and the success
rates of every executed system. Both can be given in the same call.
`--git-url <repo>` clones the set from a git repository (GitHub or any git URL)
and runs the same workflow on the clone. `--run-dir <dir>` runs TOMLs that already
exist (e.g. the `toml/` folder of a dataset repository) instead of generated ones.

Examples (from the root of a dockthor_coreset clone; DOCKTHOR is the dockthor-suite checkout):
    python3 scripts/qualification_set.py --run-dir toml -o steadystate --run-path $DOCKTHOR/build/bin/dockthor \\
        --mpi --np 24 --dockthor-args='--ga-nrun 24 --success-rate-RMSD 3'
    python3 scripts/qualification_set.py ~/vcoreset-fix --generate -o ~/results/vcoreset
    python3 scripts/qualification_set.py ~/vcoreset-fix --run-path build/bin/dockthor --mpi --np 8
    python3 scripts/qualification_set.py ~/vcoreset-fix --generate -o out --set ga-nrun=30 \\
        --run-path build/bin/dockthor --systems 1bcu 10gs
    python3 scripts/qualification_set.py ~/vcoreset-fix --generate -o steadystate --output-toml-path tomls/ss \\
        --set success-rate-RMSD=2.0 --run-path build/bin/dockthor
    python3 scripts/qualification_set.py --git-url https://github.com/<org>/<dataset>.git --generate -o steadystate \\
        --set success-rate-RMSD=2.0 --run-path build/bin/dockthor --mpi --np 2
    python3 scripts/qualification_set.py --git-url git@github.com:<org>/<dataset>.git --run-dir toml -o steadystate \\
        --run-path build/bin/dockthor --dockthor-args='--ga-nrun 24 --success-rate-RMSD 2.0'
    python3 scripts/qualification_set.py --run-dir ~/dockthor_coreset/toml -o steadystate --run-path build/bin/dockthor
"""

import argparse
import csv
import json
import re
import shlex
import subprocess
import sys
import time
import tomllib
from pathlib import Path


def read_floats(path: Path, count: int) -> list[float]:
    values = [float(v) for v in path.read_text().split()]
    if len(values) != count:
        raise ValueError(f"{path}: expected {count} value(s), found {len(values)}")
    return values


def toml_value(value) -> str:
    """Format a Python value as a TOML literal."""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return repr(value)
    if isinstance(value, list):
        return "[" + ",".join(toml_value(v) for v in value) + "]"
    # JSON string escaping is valid for TOML basic strings.
    return json.dumps(str(value))


def parse_set_option(text: str) -> tuple[str, object]:
    """Parse KEY=VALUE; VALUE is read as a TOML literal, else kept as a string."""
    key, sep, raw = text.partition("=")
    if not sep or not key.strip():
        raise argparse.ArgumentTypeError(f"expected KEY=VALUE, got '{text}'")
    try:
        value = tomllib.loads(f"v = {raw}")["v"]
    except tomllib.TOMLDecodeError:
        value = raw
    return key.strip().lstrip("-"), value


def git(*command: str, cwd: Path | None = None) -> subprocess.CompletedProcess:
    print("    $ git " + shlex.join(command), flush=True)
    return subprocess.run(["git", "-c", "advice.detachedHead=false", *command], cwd=cwd)


def repo_name(url: str) -> str:
    """https://github.com/org/set.git, git@github.com:org/set.git or /path/set -> set"""
    return re.split(r"[/:]", url.rstrip("/"))[-1].removesuffix(".git")


def clone_dataset(args) -> Path:
    """Clone (or reuse) the repository of --git-url and return the set root inside it."""
    clone_dir = (args.clone_dir or Path(repo_name(args.git_url))).expanduser().resolve()
    if (clone_dir / ".git").is_dir():
        print(f"Using the existing clone {clone_dir}")
        if args.git_update:
            ref = args.git_ref or "HEAD"
            if git("fetch", "--depth", "1", "origin", ref, cwd=clone_dir).returncode != 0 or \
                    git("checkout", "--detach", "FETCH_HEAD", cwd=clone_dir).returncode != 0:
                sys.exit(f"[error] could not update {clone_dir} to {ref}")
    elif clone_dir.exists() and any(clone_dir.iterdir()):
        sys.exit(f"[error] {clone_dir} exists and is not a git clone; choose another --clone-dir")
    else:
        print(f"Cloning {args.git_url} into {clone_dir}")
        branch = ["--branch", args.git_ref] if args.git_ref else []
        if git("clone", "--depth", "1", *branch, args.git_url, str(clone_dir)).returncode != 0:
            if not args.git_ref:
                sys.exit(f"[error] could not clone {args.git_url}")
            # --branch takes only branches and tags: a commit needs a full clone and a checkout.
            if clone_dir.exists() and any(clone_dir.iterdir()):
                sys.exit(f"[error] could not clone {args.git_url}; {clone_dir} was left behind, remove it")
            if git("clone", args.git_url, str(clone_dir)).returncode != 0 or \
                    git("checkout", "--detach", args.git_ref, cwd=clone_dir).returncode != 0:
                sys.exit(f"[error] could not clone {args.git_url} at {args.git_ref}")
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=clone_dir, capture_output=True, text=True).stdout.strip()
    print(f"Dataset commit: {commit}")
    return clone_dir / args.git_subdir if args.git_subdir else clone_dir


def infer_dataset(toml_dir: Path) -> Path | None:
    """The set root of a TOML folder: the first of toml_dir and its parents where the ligand and
    receptor of a TOML exist (their paths are relative to the set root)."""
    for toml_path in sorted(toml_dir.glob("*.toml")):
        try:
            config = tomllib.loads(toml_path.read_text())
        except (OSError, tomllib.TOMLDecodeError):
            continue
        if "ligand" not in config or "receptor" not in config:
            continue
        for root in (toml_dir, *toml_dir.parents):
            if (root / config["ligand"]).exists() and (root / config["receptor"]).exists():
                return root
        return None
    return None


def discover_systems(root: Path) -> list[str]:
    return sorted(p.name.removesuffix(".grid_center.txt") for p in (root / "config").glob("*.grid_center.txt"))


def build_config(root: Path, pdb: str, output_dir: Path, use_cofactors: bool) -> dict:
    """Collect the DockThor options of one system. Input paths are relative to the set root."""
    config_dir = root / "config"
    ligand = Path("ligand") / f"{pdb}_ligand.top"
    receptor = Path("receptor") / f"{pdb}.in"
    for path in (ligand, receptor):
        if not (root / path).is_file():
            raise FileNotFoundError(f"missing {root / path}")

    config = {
        "output-dir": str(output_dir),
        "ligand": str(ligand),
        "receptor": str(receptor),
    }
    cofactors = Path("cofactor") / pdb
    if use_cofactors and any((root / cofactors).glob("*.top")):
        config["cofactors-path"] = str(cofactors)
    config["grid-size"] = read_floats(config_dir / f"{pdb}.grid_size.txt", 3)
    config["grid-center"] = read_floats(config_dir / f"{pdb}.grid_center.txt", 3)
    config["grid-rstep"] = read_floats(config_dir / f"{pdb}.rstep.txt", 1)[0]
    return config


def write_toml(path: Path, pdb: str, config: dict) -> None:
    lines = [f"## DockThor qualification set: {pdb}", ""]
    lines += [f"{key}={toml_value(value)}" for key, value in config.items()]
    path.write_text("\n".join(lines) + "\n")


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


# The repository this script lives in (scripts/..): the default dataset when none is given.
REPOSITORY_ROOT = Path(__file__).resolve().parent.parent

ANSI_ESCAPE = re.compile(r"\x1b\[[0-9;]*m")


def system_dir_of(toml_path: Path, root: Path, output_dir: Path | None) -> Path | None:
    """DockThor writes the results of a receptor to <output-dir>/<receptor file stem>.

    output_dir overrides the TOML's output-dir (it is passed to DockThor as --output-dir)."""
    try:
        config = tomllib.loads(toml_path.read_text())
    except (OSError, tomllib.TOMLDecodeError):
        return None
    out = output_dir if output_dir is not None else config.get("output-dir")
    if out is None or "receptor" not in config:
        return None
    return root / out / Path(config["receptor"]).stem


def run_logged(command: list[str], cwd: Path, log_path: Path, echo: bool) -> int:
    """Run command, saving its console output (stdout + stderr, without colors) to log_path."""
    with open(log_path, "w") as log, subprocess.Popen(
        command, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors="replace"
    ) as process:
        for line in process.stdout:
            log.write(ANSI_ESCAPE.sub("", line))
            if echo:
                sys.stdout.write(line)
                sys.stdout.flush()
    return process.returncode


def read_success_rate(out_dir: Path) -> dict[float, tuple[float, float]]:
    """{cutoff: (SR_energy, SR_RMSD)} from DockThor's <ligand>_docked_success_rate.csv; empty if absent."""
    table = {}
    for path in out_dir.glob("*_docked_success_rate.csv"):
        with open(path, newline="") as stream:
            for row in csv.reader(stream):
                if not row or row[0].startswith("#"):
                    continue
                table[float(row[0])] = (float(row[1]), float(row[2]))
    return table


def cutoff_label(cutoff: float) -> str:
    return f"{cutoff:g}"  # 0, 0.5, 1, 1.5, ...


def write_summary(path: Path, results: list[dict]) -> None:
    """One row per executed system: complex, execution_time, RMSD_<c>..., Energy_<c>...

    Cutoffs are the union over all systems; a rate missing for a system is left empty."""
    cutoffs = sorted({c for result in results for c in result["success"]})
    header = ["complex", "execution_time"]
    header += [f"RMSD_{cutoff_label(c)}" for c in cutoffs]
    header += [f"Energy_{cutoff_label(c)}" for c in cutoffs]
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        for result in results:
            success = result["success"]
            row = [result["complex"], f"{result['time']:.2f}"]
            row += [f"{success[c][1]:.2f}" if c in success else "" for c in cutoffs]
            row += [f"{success[c][0]:.2f}" if c in success else "" for c in cutoffs]
            writer.writerow(row)


def summary_path_of(args, out_dirs: list[Path], toml_dir: Path) -> Path:
    """--summary-csv, else <common output-dir>/qualification_summary.csv, else next to the TOMLs."""
    if args.summary_csv is not None:
        return args.summary_csv.expanduser().resolve()
    parents = {d.parent for d in out_dirs}
    base = parents.pop() if len(parents) == 1 else toml_dir
    return base / "qualification_summary.csv"


def run_systems(args) -> list[str]:
    """--systems, else every <pdb>.toml of --run-dir, else every system of config/."""
    if args.systems:
        return args.systems
    if args.run_toml_dir is not None:
        return sorted(p.stem for p in args.run_toml_dir.glob("*.toml"))
    return discover_systems(args.dataset)


def run(args, systems: list[str]) -> int:
    dockthor = args.run_path.expanduser().resolve()
    if not dockthor.is_file():
        print(f"[error] DockThor binary not found: {dockthor}", file=sys.stderr)
        return 1

    launcher = []
    if args.mpi:
        launcher = [args.mpirun, "-n", str(args.np), *shlex.split(args.mpi_args)]
    dockthor_args = shlex.split(args.dockthor_args)
    output_dir = args.output_dir.expanduser() if args.output_dir is not None else None
    if output_dir is not None:
        dockthor_args = ["--output-dir", str(output_dir), *dockthor_args]

    toml_dir = args.run_toml_dir or args.toml_dir
    tomls = [toml_dir / f"{pdb}.toml" for pdb in systems]
    missing = [t.name for t in tomls if not t.is_file()]
    if missing:
        print(f"[warn] {len(missing)} TOML file(s) not found in {toml_dir}: {' '.join(missing[:10])}")
    if not systems:
        print(f"[warn] no TOML file to run in {toml_dir}")
    tomls = [t for t in tomls if t.is_file()]

    failed = []
    results = []
    out_dirs = [d for t in tomls if (d := system_dir_of(t, args.dataset, output_dir)) is not None]
    summary_path = summary_path_of(args, out_dirs, toml_dir)
    total_start = time.monotonic()
    for index, toml_path in enumerate(tomls, start=1):
        pdb = toml_path.stem
        out_dir = system_dir_of(toml_path, args.dataset, output_dir)
        if args.skip_done and out_dir is not None and any(out_dir.glob("*_docked.csv")):
            print(f"[{index}/{len(tomls)}] {pdb}: skip, results already in {out_dir}")
            continue

        # Run from the set root: the input paths in the TOML are relative to it.
        command = [*launcher, str(dockthor), "--config", str(toml_path), *dockthor_args]
        print(f"[{index}/{len(tomls)}] {pdb}: {shlex.join(command)}", flush=True)
        if args.dry_run:
            continue

        start = time.monotonic()
        if out_dir is None:
            print("    [warn] no output-dir/receptor in the TOML, console output not saved")
            log_path = None
            returncode = subprocess.run(command, cwd=args.dataset).returncode
        else:
            out_dir.mkdir(parents=True, exist_ok=True)
            log_path = out_dir / f"{out_dir.name}_run.log"
            returncode = run_logged(command, args.dataset, log_path, args.verbose)
        elapsed = time.monotonic() - start
        status = "ok" if returncode == 0 else f"FAILED (exit {returncode})"
        print(f"    {status} in {elapsed:.1f}s" + (f", log: {log_path}" if log_path else ""))
        results.append({"complex": pdb,
                        "time": elapsed,
                        "success": read_success_rate(out_dir) if returncode == 0 and out_dir else {}})
        write_summary(summary_path, results)  # rewritten after every system, so an interrupted set keeps its rows
        if returncode != 0:
            failed.append(pdb)
            if args.stop_on_error:
                break

    print(f"Ran {len(tomls)} system(s) in {time.monotonic() - total_start:.1f}s, {len(failed)} failed.")
    if results:
        print(f"Summary: {summary_path}")
        if not any(r["success"] for r in results):
            print("[warn] no success rates found; set success-rate-RMSD (--set or --dockthor-args) to get them")
    if failed:
        print("Failed: " + " ".join(failed))
    return len(failed)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate and run DockThor TOML configurations for a qualification set.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("dataset", type=Path, nargs="?",
                        help="Qualification set root (contains config/, ligand/, receptor/, ...). Omit it with "
                             "--git-url, or with --run-dir to infer it from the TOMLs; with none of them, the "
                             "repository holding this script is used.")
    parser.add_argument("--systems", nargs="+", metavar="PDB", help="Only these PDB IDs (default: all in config/).")

    src = parser.add_argument_group("dataset from git")
    src.add_argument("--git-url", metavar="URL",
                     help="Clone the qualification set from this git repository (e.g. "
                          "https://github.com/<org>/<repo>.git) and run the workflow on the clone.")
    src.add_argument("--git-ref", metavar="REF", help="Branch, tag or commit to check out (default: the default branch).")
    src.add_argument("--git-subdir", metavar="PATH",
                     help="Set root inside the repository, if it is not the repository root.")
    src.add_argument("--clone-dir", type=Path, metavar="DIR",
                     help="Where to clone (default: ./<repo name>). An existing clone there is reused as is.")
    src.add_argument("--git-update", action="store_true",
                     help="Update an existing clone to --git-ref (or the remote default branch) before running.")

    gen = parser.add_argument_group("generation")
    gen.add_argument("--generate", action="store_true", help="Write the <pdb>.toml files.")
    gen.add_argument("--output-toml-path", type=Path, metavar="DIR",
                     help="Directory of the TOML files, written by --generate and read by --run-path "
                          "(default: the dataset root). Relative to the current directory. The ligand, receptor "
                          "and cofactor paths inside the TOMLs stay relative to the dataset root.")
    gen.add_argument("-o", "--output-dir", type=Path,
                     help="DockThor output directory; DockThor writes each system to <output-dir>/<pdb>. "
                          "Written into the TOMLs with --generate; with --run-path it is also passed to DockThor "
                          "as --output-dir, overriding the TOMLs. A relative path is relative to the dataset root.")
    gen.add_argument("--set", action="append", default=[], type=parse_set_option, metavar="KEY=VALUE",
                     help="Extra DockThor option written to every TOML, e.g. --set ga-nrun=30 "
                          "--set seed=7 --set 'ga-evolution-type=\"GENERATIONAL_OCL\"'. Repeatable.")
    gen.add_argument("--no-cofactors", action="store_true", help="Do not add cofactors-path.")
    gen.add_argument("--overwrite", action="store_true", help="Overwrite existing TOML files.")

    exe = parser.add_argument_group("execution")
    exe.add_argument("--run-path", type=Path, metavar="DOCKTHOR", help="Path to the dockthor binary; runs every TOML.")
    exe.add_argument("--run-dir", type=Path, metavar="DIR",
                     help="Run the TOMLs that already exist in DIR (every <pdb>.toml, or --systems) instead of "
                          "the generated ones. Their input paths must be relative to the dataset root. With a "
                          "dataset or --git-url, a relative DIR is relative to the dataset root (e.g. 'toml'); "
                          "alone, DIR is relative to the current directory and the dataset root is inferred "
                          "from the TOMLs (DIR or the first parent holding their ligand and receptor).")
    exe.add_argument("--mpi", action="store_true", default=False, help="Run DockThor with MPI.")
    exe.add_argument("--np", type=int, default=4, help="Number of MPI processes (with --mpi).")
    exe.add_argument("--mpirun", default="mpirun", help="MPI launcher command.")
    exe.add_argument("--mpi-args", default="", help="Extra launcher arguments; use the = form, e.g. --mpi-args=--oversubscribe.")
    exe.add_argument("--dockthor-args", default="", help="Extra command-line arguments passed to DockThor; use the = form, e.g. --dockthor-args='--seed 7'.")
    exe.add_argument("--summary-csv", type=Path, metavar="FILE",
                     help="Summary of the executed systems: complex, execution_time (s), RMSD_<cutoff>... and "
                          "Energy_<cutoff>... success rates in %% (default: qualification_summary.csv in the common "
                          "output directory, else next to the TOMLs).")
    exe.add_argument("--skip-done", action="store_true", help="Skip systems that already have a *_docked.csv in <output-dir>/<pdb>.")
    exe.add_argument("--stop-on-error", action="store_true", help="Stop at the first failing system.")
    exe.add_argument("--dry-run", action="store_true", help="Print the commands without running them.")
    exe.add_argument("-v", "--verbose", action="store_true",
                     help="Also print DockThor output to the terminal; it is always saved to "
                          "<output-dir>/<pdb>/<pdb>_run.log.")

    args = parser.parse_args(argv)
    if args.dataset is not None and args.git_url is not None:
        parser.error("give either the dataset directory or --git-url, not both")
    if args.dataset is None and args.git_url is None and args.run_dir is None:
        if not (REPOSITORY_ROOT / "config").is_dir():
            parser.error("give the dataset directory, --git-url, or --run-dir with the TOMLs")
        args.dataset = REPOSITORY_ROOT  # the script runs from a dataset repository: use it
    if args.dataset is None and args.git_url is None and args.generate and args.run_dir is not None:
        parser.error("--generate needs the dataset directory or --git-url")
    if args.git_url is None and (args.git_ref or args.git_subdir or args.clone_dir or args.git_update):
        parser.error("--git-ref, --git-subdir, --clone-dir and --git-update need --git-url")
    if not args.generate and args.run_path is None and args.git_url is None:
        parser.error("nothing to do: pass --generate and/or --run-path")
    if args.generate and args.output_dir is None:
        parser.error("--generate requires --output-dir")
    if args.run_dir is not None and args.run_path is None:
        parser.error("--run-dir requires --run-path")
    args.parser = parser
    return args


def resolve_dataset(args) -> None:
    """Set args.dataset (cloning it first with --git-url, or inferring it from --run-dir) and args.toml_dir."""
    if args.dataset is None and args.git_url is None:
        toml_dir = args.run_dir.expanduser().resolve()
        if not toml_dir.is_dir():
            args.parser.error(f"--run-dir {toml_dir} is not a directory")
        dataset = infer_dataset(toml_dir)
        if dataset is None:
            args.parser.error(f"cannot find the dataset root of the TOMLs in {toml_dir}: no folder from there up "
                              "holds their ligand and receptor; pass the dataset directory")
        print(f"Dataset root (from the TOMLs): {dataset}")
        args.run_dir = toml_dir  # absolute, so it does not depend on the dataset root below
    else:
        dataset = clone_dataset(args) if args.git_url else args.dataset
    args.dataset = dataset.expanduser().resolve()
    args.toml_dir = args.output_toml_path.expanduser().resolve() if args.output_toml_path else args.dataset
    args.run_toml_dir = None
    if args.run_dir is not None:
        args.run_toml_dir = (args.dataset / args.run_dir.expanduser()).resolve()  # an absolute DIR stays as is
        if not args.run_toml_dir.is_dir():
            args.parser.error(f"--run-dir {args.run_toml_dir} is not a directory")
    needs_config = args.generate or (args.run_path is not None and args.run_dir is None and not args.systems)
    if needs_config and not (args.dataset / "config").is_dir():
        args.parser.error(f"{args.dataset} has no config/ directory")


def main(argv=None) -> int:
    args = parse_args(argv)
    resolve_dataset(args)
    errors = 0
    if args.generate:
        errors += generate(args, args.systems or discover_systems(args.dataset))
    if args.run_path is not None:
        errors += run(args, run_systems(args))
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
