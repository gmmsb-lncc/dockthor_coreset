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
overall_min_RMSD <= cutoff, column overall) and every cutoff 0.5, 1.0, ... up to `--rmsd-cutoff`. Both can be given in the same call.
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
    python3 scripts/qualification_set.py --git-url git@github.com:<org>/<dataset>.git --config-dir toml -o steadystate \\
        --dockthor-path build/bin/dockthor --dockthor-args='--ga-nrun 24 --success-rate --success-rate-by-run 2.0'
    python3 scripts/qualification_set.py --config-dir ~/dockthor_coreset/toml -o steadystate --dockthor-path build/bin/dockthor \\
        --rmsd-cutoff 3
    python3 scripts/qualification_set.py --generate-summary steadystate/qualification_summary.csv --rmsd-cutoff 3
    python3 scripts/qualification_set.py --collect ci_collect/steadystate_symmetric --config-dir toml
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


SUCCESS_RATE_CSV = "_docked_success_rate.csv"  # suffix of DockThor's --success-rate report, <ligand>_docked_success_rate.csv


def read_success_rate_file(path: Path) -> dict[str, str]:
    """{column: value} of one <ligand>_docked_success_rate.csv; empty if it is malformed."""
    with open(path, newline="") as stream:
        rows = [row for row in csv.reader(stream) if row]
    return dict(zip(rows[0], rows[1])) if len(rows) >= 2 else {}


def read_success_rate(out_dir: Path) -> dict[str, str]:
    """{column: value} of DockThor's one-row <ligand>_docked_success_rate.csv (--success-rate):
    complex, execution_time(s), num_top, overall_min_RMSD_E, overall_min_RMSD, top1_E, top1_RMSD, ..., topN_E,
    topN_RMSD. Any N is accepted (overall_min_RMSD* is absent in files of older DockThor builds); empty if absent
    or malformed."""
    for path in sorted(out_dir.glob(f"*{SUCCESS_RATE_CSV}")):
        if report := read_success_rate_file(path):
            return report
    return {}


def column_key(column: str) -> tuple:
    """Summary column order: complex, execution_time(s), num_top, overall_min_RMSD_E, overall_min_RMSD, top1_E,
    top1_RMSD, top2_E, ... Unknown columns go last, by name."""
    fixed = ["complex", "execution_time(s)", "num_top", "overall_min_RMSD_E", "overall_min_RMSD"]
    if column in fixed:
        return (0, fixed.index(column))
    if match := re.fullmatch(r"top(\d+)_(E|RMSD)", column):
        return (1, int(match[1]), match[2] != "E")
    return (2, column)


def top_ks(reports: list[dict]) -> list[int]:
    """The K of every topK_RMSD column in the reports, ascending."""
    return sorted({int(m[1]) for report in reports for c in report if (m := re.fullmatch(r"top(\d+)_RMSD", c))})


def write_summary(path: Path, reports: list[dict]) -> None:
    """qualification_summary.csv: the rows of every <ligand>_docked_success_rate.csv, one per complex.

    Columns are the union over all complexes, so complexes with different N fit in one file; a value
    missing for a complex is left empty."""
    columns = sorted({c for report in reports for c in report}, key=column_key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        for report in reports:
            writer.writerow([report.get(c, "") for c in columns])


def success_rate_cutoffs(max_RMSD: float) -> list[float]:
    """0.5, 1.0, ... up to max_RMSD rounded up to the 0.5 grid (0.3 -> [0.5], 2.2 -> [..., 2.5]),
    the cutoffs of DockThor's --success-rate-by-run."""
    cutoffs = []
    while not cutoffs or cutoffs[-1] < max_RMSD - 1e-4:
        cutoffs.append(0.5 * (len(cutoffs) + 1))
    return cutoffs


OVERALL = "overall"  # success-rate column of overall_min_RMSD, before the topK ones


def rate_columns(reports: list[dict]) -> list[tuple[str, str]]:
    """(label, summary column) of every success_rate.txt column: overall (overall_min_RMSD, the lowest RMSD over
    all the leaders) when a report has it, then top1 ... topN."""
    columns = [(OVERALL, "overall_min_RMSD")] if any(report.get("overall_min_RMSD") for report in reports) else []
    return columns + [(f"top{k}", f"top{k}_RMSD") for k in top_ks(reports)]


def success_rates(reports: list[dict], cutoffs: list[float]) -> dict[tuple[float, str], tuple[int, int]]:
    """{(cutoff, label): (complexes with that RMSD <= cutoff, complexes with that RMSD)}, label topK or overall."""
    table = {}
    for label, column in rate_columns(reports):
        values = [float(report[column]) for report in reports if report.get(column)]
        for cutoff in cutoffs:
            table[cutoff, label] = (sum(value <= cutoff for value in values), len(values))
    return table


def write_success_rate_txt(path: Path, reports: list[dict], cutoffs: list[float], num_run: int) -> None:
    """success_rate.txt: for every cutoff c and every topK, the % of complexes whose topK_RMSD <= c, and in the
    overall column the % whose overall_min_RMSD <= c, over the complexes with a _success_rate.csv."""
    columns = rate_columns(reports)
    table = success_rates(reports, cutoffs)
    lines = ["Success rate (%) of the top-N report (DockThor --success-rate), topK_RMSD <= cutoff"]
    if any(label == OVERALL for label, _ in columns):
        lines.append("overall: overall_min_RMSD <= cutoff (the lowest-RMSD leader of all the leaders)")
    lines += [f"Complexes with results: {len(reports)} of {num_run}", ""]
    lines.append(f"{'RMSD_cutoff':>11}" + "".join(f"{label:>9}" for label, _ in columns))
    for cutoff in cutoffs:
        cells = []
        for label, _ in columns:
            hits, total = table[cutoff, label]
            cells.append(f"{100 * hits / total:9.2f}" if total else f"{'N/A':>9}")
        lines.append(f"{cutoff:>11g}" + "".join(cells))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def read_summary(path: Path) -> list[dict]:
    """The rows of a qualification_summary.csv, one {column: value} per complex."""
    with open(path, newline="") as stream:
        return [row for row in csv.DictReader(stream) if any(row.values())]


def generate_summary(args) -> int:
    """--generate-summary: success_rate.txt, in the current directory, from an existing summary CSV."""
    summary_path = args.generate_summary.expanduser().resolve()
    try:
        reports = read_summary(summary_path)
    except OSError as error:
        print(f"[error] {error}", file=sys.stderr)
        return 1
    if not top_ks(reports):
        print(f"[error] {summary_path}: no topK_RMSD column, not a qualification_summary.csv", file=sys.stderr)
        return 1
    success_rate_path = Path.cwd() / "success_rate.txt"
    write_success_rate_txt(success_rate_path, reports, success_rate_cutoffs(args.rmsd_cutoff), len(reports))
    print(f"Success rates of {len(reports)} complex(es) from {summary_path}: {success_rate_path}")
    return 0


def collect(args) -> int:
    """--collect DIR: qualification_summary.csv and success_rate.txt from every <ligand>_docked_success_rate.csv
    under DIR (flat, or in <pdb>/ folders), e.g. the per-complex reports of parallel runs gathered into one folder.
    With --config-dir, the complexes are counted against its TOMLs and the missing ones are listed."""
    root = args.collect.expanduser().resolve()
    if not root.is_dir():
        print(f"[error] --collect {root} is not a directory", file=sys.stderr)
        return 1

    reports = []
    seen = set()
    for path in sorted(root.rglob(f"*{SUCCESS_RATE_CSV}")):
        report = read_success_rate_file(path)
        if not report:
            print(f"[warn] {path}: malformed, skipped")
            continue
        name = report.get("complex") or path.name.removesuffix(SUCCESS_RATE_CSV)
        if name in seen:
            print(f"[warn] {path}: a second report of {name}, skipped")
            continue
        seen.add(name)
        reports.append(report)
    reports.sort(key=lambda report: report.get("complex", ""))
    if not reports or not top_ks(reports):
        print(f"[error] no *{SUCCESS_RATE_CSV} with topK_RMSD columns under {root}", file=sys.stderr)
        return 1

    expected = []
    if args.config_dir is not None:
        toml_dir = args.config_dir.expanduser().resolve()
        if not toml_dir.is_dir():
            print(f"[error] --config-dir {toml_dir} is not a directory", file=sys.stderr)
            return 1
        expected = sorted(p.stem for p in toml_dir.glob("*.toml"))
    num_run = len(expected) if expected else len(reports)

    summary_path = args.summary_csv.expanduser().resolve() if args.summary_csv else root / "qualification_summary.csv"
    success_rate_path = summary_path.with_name("success_rate.txt")
    write_summary(summary_path, reports)
    write_success_rate_txt(success_rate_path, reports, success_rate_cutoffs(args.rmsd_cutoff), num_run)
    print(f"Collected {len(reports)} of {num_run} complex(es) from {root}")
    print(f"Summary: {summary_path}")
    print(f"Success rates: {success_rate_path}")
    if missing := [pdb for pdb in expected if pdb not in seen]:
        print(f"[warn] {len(missing)} complex(es) without results: {' '.join(missing)}")
    return 0


def summary_path_of(args, out_dirs: list[Path], toml_dir: Path) -> Path:
    """--summary-csv, else <common output-dir>/qualification_summary.csv, else next to the TOMLs."""
    if args.summary_csv is not None:
        return args.summary_csv.expanduser().resolve()
    parents = {d.parent for d in out_dirs}
    base = parents.pop() if len(parents) == 1 else toml_dir
    return base / "qualification_summary.csv"


def run_systems(args) -> list[str]:
    """--systems, else every <pdb>.toml of --config-dir, else every system of config/."""
    if args.systems:
        return args.systems
    if args.run_toml_dir is not None:
        return sorted(p.stem for p in args.run_toml_dir.glob("*.toml"))
    return discover_systems(args.dataset)


def run(args, systems: list[str]) -> int:
    dockthor = args.dockthor_path.expanduser().resolve()
    if not dockthor.is_file():
        print(f"[error] DockThor binary not found: {dockthor}", file=sys.stderr)
        return 1

    launcher = []
    if args.mpi:
        launcher = [args.mpirun, "-n", str(args.np), *shlex.split(args.mpi_args)]
    dockthor_args = shlex.split(args.dockthor_args)
    if "--success-rate" not in dockthor_args:
        # Without it DockThor writes no _success_rate.csv, and the docking is lost for the summary.
        dockthor_args.append("--success-rate")
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
    reports = []  # one _success_rate.csv row per complex with results, run now or skipped as done
    out_dirs = [d for t in tomls if (d := system_dir_of(t, args.dataset, output_dir)) is not None]
    summary_path = summary_path_of(args, out_dirs, toml_dir)
    success_rate_path = summary_path.with_name("success_rate.txt")
    cutoffs = success_rate_cutoffs(args.rmsd_cutoff)

    def write_reports() -> None:  # rewritten after every complex, so an interrupted set keeps what is done
        write_summary(summary_path, reports)
        write_success_rate_txt(success_rate_path, reports, cutoffs, len(tomls))

    total_start = time.monotonic()
    for index, toml_path in enumerate(tomls, start=1):
        pdb = toml_path.stem
        out_dir = system_dir_of(toml_path, args.dataset, output_dir)
        if args.skip_done and out_dir is not None and any(out_dir.glob("*_docked.csv")):
            print(f"[{index}/{len(tomls)}] {pdb}: skip, results already in {out_dir}")
            if report := read_success_rate(out_dir):
                reports.append(report)
                if not args.dry_run:
                    write_reports()
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
        if returncode == 0 and out_dir is not None and (report := read_success_rate(out_dir)):
            reports.append(report)
        write_reports()
        if returncode != 0:
            failed.append(pdb)
            if args.stop_on_error:
                break

    print(f"Ran {len(tomls)} system(s) in {time.monotonic() - total_start:.1f}s, {len(failed)} failed.")
    if not args.dry_run and tomls:
        print(f"Summary: {summary_path} ({len(reports)} complex(es) with a _success_rate.csv)")
        print(f"Success rates: {success_rate_path}")
        if not reports:
            print("[warn] no _success_rate.csv found; check the run logs")
    if failed:
        print("Failed: " + " ".join(failed))
    return len(failed)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate and run DockThor TOML configurations for a qualification set.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("dataset", type=Path, nargs="?",
                        help="Qualification set root (contains config/, ligand/, receptor/, ...). "
                             "Omit it with --git-url, or with --config-dir to infer it from the TOMLs.")
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
                     help="Directory of the TOML files, written by --generate and read by --dockthor-path "
                          "(default: the dataset root). Relative to the current directory. The ligand, receptor "
                          "and cofactor paths inside the TOMLs stay relative to the dataset root.")
    gen.add_argument("-o", "--output-dir", type=Path,
                     help="DockThor output directory; DockThor writes each system to <output-dir>/<pdb>. "
                          "Written into the TOMLs with --generate; with --dockthor-path it is also passed to DockThor "
                          "as --output-dir, overriding the TOMLs. A relative path is relative to the dataset root.")
    gen.add_argument("--set", action="append", default=[], type=parse_set_option, metavar="KEY=VALUE",
                     help="Extra DockThor option written to every TOML, e.g. --set ga-nrun=30 "
                          "--set seed=7 --set 'ga-evolution-type=\"GENERATIONAL_OCL\"'. Repeatable.")
    gen.add_argument("--no-cofactors", action="store_true", help="Do not add cofactors-path.")
    gen.add_argument("--overwrite", action="store_true", help="Overwrite existing TOML files.")

    exe = parser.add_argument_group("execution")
    exe.add_argument("--dockthor-path", type=Path, metavar="DOCKTHOR", help="Path to the dockthor binary; runs every TOML.")
    exe.add_argument("--config-dir", type=Path, metavar="DIR",
                     help="Run the TOMLs that already exist in DIR (every <pdb>.toml, or --systems) instead of "
                          "the generated ones. Their input paths must be relative to the dataset root. With a "
                          "dataset or --git-url, a relative DIR is relative to the dataset root (e.g. 'toml'); "
                          "alone, DIR is relative to the current directory and the dataset root is inferred "
                          "from the TOMLs (DIR or the first parent holding their ligand and receptor).")
    exe.add_argument("--mpi", action="store_true", default=False, help="Run DockThor with MPI.")
    exe.add_argument("--np", type=int, default=4, help="Number of MPI processes (with --mpi).")
    exe.add_argument("--mpirun", default="mpirun", help="MPI launcher command.")
    exe.add_argument("--mpi-args", default="", help="Extra launcher arguments; use the = form, e.g. --mpi-args=--oversubscribe.")
    exe.add_argument("--dockthor-args", default="", help="Extra command-line arguments passed to DockThor; use the = form, e.g. --dockthor-args='--seed 7'. "
                          "--success-rate is always added.")
    exe.add_argument("--summary-csv", type=Path, metavar="FILE",
                     help="Consolidation of every <pdb>/<ligand>_docked_success_rate.csv (DockThor --success-rate): "
                          "complex, execution_time(s), num_top, overall_min_RMSD_E, overall_min_RMSD (energy and RMSD "
                          "of the lowest-RMSD leader of the all-runs clustering), top1_E, top1_RMSD, ..., "
                          "topN_E, topN_RMSD, one row "
                          "per complex (default: qualification_summary.csv in the common output directory, else "
                          "next to the TOMLs). success_rate.txt is written next to it.")
    exe.add_argument("--rmsd-cutoff", type=float, default=2.0, metavar="RMSD",
                     help="Largest RMSD cutoff of success_rate.txt: for each cutoff 0.5, 1.0, ... up to RMSD "
                          "(rounded up to 0.5) and each topK, the %% of complexes with topK_RMSD <= cutoff, and in the "
                          "overall column with overall_min_RMSD <= cutoff.")
    exe.add_argument("--generate-summary", type=Path, metavar="CSV",
                     help="Only write success_rate.txt, in the current directory, from the existing "
                          "qualification_summary.csv CSV (with --rmsd-cutoff); nothing is generated or docked.")
    exe.add_argument("--collect", type=Path, metavar="DIR",
                     help="Only consolidate: write qualification_summary.csv (or --summary-csv) and success_rate.txt "
                          "from every <ligand>_docked_success_rate.csv under DIR, flat or in <pdb>/ folders (e.g. "
                          "the reports of parallel runs gathered in one folder). With --config-dir, the complexes "
                          "are counted against its TOMLs and the missing ones listed. Nothing is generated or docked.")
    exe.add_argument("--skip-done", action="store_true", help="Skip systems that already have a *_docked.csv in <output-dir>/<pdb>.")
    exe.add_argument("--stop-on-error", action="store_true", help="Stop at the first failing system.")
    exe.add_argument("--dry-run", action="store_true", help="Print the commands without running them.")
    exe.add_argument("-v", "--verbose", action="store_true",
                     help="Also print DockThor output to the terminal; it is always saved to "
                          "<output-dir>/<pdb>/<pdb>_run.log.")

    args = parser.parse_args(argv)
    if args.collect is not None:
        if (args.dataset is not None or args.git_url or args.generate or args.dockthor_path
                or args.generate_summary is not None):
            parser.error("--collect runs alone (with --config-dir, --summary-csv and --rmsd-cutoff): drop the "
                         "dataset, --git-url, --generate, --dockthor-path and --generate-summary")
        if args.rmsd_cutoff <= 0:
            parser.error("--rmsd-cutoff must be positive")
        args.parser = parser
        return args
    if args.generate_summary is not None:
        if args.dataset is not None or args.git_url or args.config_dir or args.generate or args.dockthor_path:
            parser.error("--generate-summary runs alone: drop the dataset, --git-url, --config-dir, "
                         "--generate and --dockthor-path")
        if args.rmsd_cutoff <= 0:
            parser.error("--rmsd-cutoff must be positive")
        args.parser = parser
        return args
    if args.dataset is not None and args.git_url is not None:
        parser.error("give either the dataset directory or --git-url, not both")
    if args.dataset is None and args.git_url is None and args.config_dir is None:
        parser.error("give the dataset directory, --git-url, or --config-dir with the TOMLs "
                     "(or --generate-summary CSV)")
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
            args.parser.error(f"cannot find the dataset root of the TOMLs in {toml_dir}: no folder from there up "
                              "holds their ligand and receptor; pass the dataset directory")
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
        return collect(args)
    if args.generate_summary is not None:
        return generate_summary(args)
    resolve_dataset(args)
    errors = 0
    if args.generate:
        errors += generate(args, args.systems or discover_systems(args.dataset))
    if args.dockthor_path is not None:
        errors += run(args, run_systems(args))
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
