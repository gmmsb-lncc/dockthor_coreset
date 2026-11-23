#!/usr/bin/env python3
"""
Rescore a finished docking run with DockTScore, without docking anything again.

`qualification_set.py --docktscore-path` rescores while it docks (README 4.7); this script does the same work on a
run that is already on disk -- a run docked before the rescoring existed, one docked with other models, or one
docked elsewhere -- and rewrites its reports:

    <docked-dir>/<pdb>/<ligand>_docked.mol2        what it reads: the poses of every complex
    <docked-dir>/<pdb>/<pdb>_top1_<model>.score    what DockTScore writes, per model
    <docked-dir>/<pdb>/<pdb>_scores.csv            one row per complex: its score per model
    <docked-dir>/qualification_summary.csv         rewritten: the docktscore_<model> columns and the affinities
    <docked-dir>/scoring_power.txt / .csv          the statistics per model (README 5.5)

The summary keeps everything DockThor reported: the rows come from the `<ligand>_docked_success_rate.csv` of each
complex, as `summarize_results.py` reads them, and only the score and affinity columns are added. A complex is
skipped (with a warning) when it has no poses; one already scored is scored again unless --skip-done is given.

The set the run docked is needed for the receptors, the cofactors and `index/affinity.csv`: give it as --dataset,
or let it be inferred from the TOMLs with --config-dir.

Examples:
    python3 scripts/rescore_results.py steadystate --dataset v2013-coreset \\
        --docktscore-path $DOCKTHOR/build/bin/docktscore
    python3 scripts/rescore_results.py steadystate --config-dir v2013-coreset/config \\
        --docktscore-path $DOCKTHOR/build/bin/docktscore --docktscore-models general_dock protease
    python3 scripts/rescore_results.py ci_collect/steadystate_symmetric --dataset v2013-coreset \\
        --docktscore-path $DOCKTHOR/build/bin/docktscore --skip-done
"""

import argparse
import shlex
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from coreset import scoring, summary  # noqa: E402
from coreset.affinity import read_affinities  # noqa: E402
from coreset.dataset import infer_dataset  # noqa: E402

DEFAULT_SUMMARY = "qualification_summary.csv"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Rescore a finished docking run with DockTScore, without docking again.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("docked", type=Path, metavar="DIR", help="A run's output directory, one <pdb>/ per complex.")
    parser.add_argument(
        "--dataset", type=Path, help="The set that was docked (receptors, cofactors, index/affinity.csv)."
    )
    parser.add_argument(
        "--config-dir", type=Path, metavar="DIR", help="Infer the set from the TOMLs of this folder instead."
    )
    parser.add_argument(
        "--complex", nargs="+", metavar="PDB", help="Only these PDB IDs (default: every complex of DIR)."
    )
    parser.add_argument(
        "--docktscore-path", type=Path, required=True, metavar="DOCKTSCORE", help="Path to the docktscore binary."
    )
    parser.add_argument(
        "--docktscore-models",
        "--docktscore-model",
        nargs="+",
        default=scoring.DEFAULT_MODELS,
        metavar="MODEL",
        help="DockTScore coefficients to score with, one summary column each (docktscore_<model>).",
    )
    parser.add_argument(
        "--docktscore-args",
        default="",
        help="Extra arguments for DockTScore; use the = form, e.g. --docktscore-args='-v'.",
    )
    parser.add_argument(
        "--summary-csv", type=Path, metavar="FILE", help=f"The summary to write (default DIR/{DEFAULT_SUMMARY})."
    )
    parser.add_argument(
        "--rmsd-cutoff",
        type=float,
        default=3.0,
        metavar="RMSD",
        help="Largest RMSD cutoff of success_rate.txt (default 3.0).",
    )
    parser.add_argument("--skip-done", action="store_true", help="Skip a complex that already has every model's score.")
    args = parser.parse_args(argv)
    if args.dataset is None and args.config_dir is None:
        parser.error("give --dataset or --config-dir")
    if args.rmsd_cutoff <= 0:
        parser.error("--rmsd-cutoff must be positive")
    args.docktscore_args = shlex.split(args.docktscore_args)
    args.docktscore_path = args.docktscore_path.expanduser().resolve()
    if not args.docktscore_path.is_file():
        parser.error(f"--docktscore-path {args.docktscore_path} is not a file")
    args.parser = parser
    return args


def resolve_dataset(args) -> None:
    """Set args.dataset from --dataset or from the TOMLs of --config-dir."""
    if args.dataset is not None:
        args.dataset = args.dataset.expanduser().resolve()
        if not args.dataset.is_dir():
            args.parser.error(f"--dataset {args.dataset} is not a directory")
        return
    toml_dir = args.config_dir.expanduser().resolve()
    if not toml_dir.is_dir():
        args.parser.error(f"--config-dir {toml_dir} is not a directory")
    dataset = infer_dataset(toml_dir)
    if dataset is None:
        args.parser.error(
            f"cannot find the dataset root of the TOMLs in {toml_dir}: no folder from there up holds their "
            "ligand and receptor; pass --dataset"
        )
    print(f"Dataset root (from the TOMLs): {dataset}")
    args.dataset = dataset


def complexes_of(docked: Path, only: list[str] | None) -> list[tuple[str, Path]]:
    """(complex, its directory) of every complex of a run's output directory that holds poses, by name."""
    found = []
    for directory in sorted(path for path in docked.iterdir() if path.is_dir()):
        if only is not None and directory.name not in only:
            continue
        if scoring.docked_poses(directory) is not None:
            found.append((directory.name, directory))
    return found


def is_done(report: dict, columns: list[str]) -> bool:
    """True when the row already holds every model's score."""
    return all(report.get(column) for column in columns)


def main(argv=None) -> int:
    args = parse_args(argv)
    resolve_dataset(args)
    docked = args.docked.expanduser().resolve()
    if not docked.is_dir():
        args.parser.error(f"{docked} is not a directory")

    reports, _ = summary.collect_reports(docked)
    if not reports:
        print(f"[error] no *{summary.SUCCESS_RATE_CSV} under {docked}", file=sys.stderr)
        return 1
    reports = summary.merge_scores(reports, docked)  # scores already on disk, so --skip-done can see them
    by_complex = {report.get("complex", ""): report for report in reports}
    columns = [scoring.column_of(model) for model in args.docktscore_models]

    found = complexes_of(docked, args.complex)
    if not found:
        print(f"[error] no complex with poses (*_docked.mol2) in {docked}", file=sys.stderr)
        return 1

    start = time.monotonic()
    scored = 0
    for index, (pdb, out_dir) in enumerate(found, start=1):
        report = by_complex.get(pdb)
        if report is None:
            print(f"[{index}/{len(found)}] {pdb}: no success-rate report, skipped")
            continue
        if args.skip_done and is_done(report, columns):
            print(f"[{index}/{len(found)}] {pdb}: skip, already scored")
            continue
        print(f"[{index}/{len(found)}] {pdb}: {' '.join(args.docktscore_models)}", flush=True)
        if scores := scoring.score_complex(args, pdb, out_dir):
            report.update(scores)
            scored += 1

    affinities = read_affinities(args.dataset)
    reports = [summary.add_affinity(report, affinities) for report in reports]
    summary_path = (args.summary_csv or docked / DEFAULT_SUMMARY).expanduser().resolve()
    summary.write_summary(summary_path, reports)
    summary.write_success_rate_txt(
        summary_path.with_name("success_rate.txt"), reports, summary.success_rate_cutoffs(args.rmsd_cutoff), len(found)
    )
    print(f"Rescored {scored} of {len(found)} complex(es) in {time.monotonic() - start:.1f}s")
    print(f"Summary: {summary_path}")
    power = summary_path.with_name("scoring_power.txt")
    if summary.write_scoring_power(power, reports):
        print(f"Scoring power: {power} and {power.with_suffix('.csv')}")
    else:
        print("[warn] no score in any row: nothing to report")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
