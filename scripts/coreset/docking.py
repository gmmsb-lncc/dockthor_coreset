"""Running DockThor on the TOMLs of the set: one complex after the other, a log per complex, and the summary and
success rates rewritten after every complex."""

import re
import shlex
import subprocess
import sys
import time
import tomllib
from pathlib import Path

from coreset.dataset import discover_systems
from coreset.summary import read_success_rate, success_rate_cutoffs, write_success_rate_txt, write_summary

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
    with (
        open(log_path, "w") as log,
        subprocess.Popen(
            command, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors="replace"
        ) as process,
    ):
        for line in process.stdout:
            log.write(ANSI_ESCAPE.sub("", line))
            if echo:
                sys.stdout.write(line)
                sys.stdout.flush()
    return process.returncode


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


def dockthor_command_args(args) -> tuple[list[str], list[str], Path | None]:
    """(launcher, DockThor arguments, --output-dir override) of the run."""
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
    return launcher, dockthor_args, output_dir


def run(args, systems: list[str]) -> int:
    """Dock every system, from the set root; the number of failed systems."""
    dockthor = args.dockthor_path.expanduser().resolve()
    if not dockthor.is_file():
        print(f"[error] DockThor binary not found: {dockthor}", file=sys.stderr)
        return 1
    launcher, dockthor_args, output_dir = dockthor_command_args(args)

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
