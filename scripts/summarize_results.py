#!/usr/bin/env python3
"""
Merge DockThor success-rate reports into qualification_summary.csv and success_rate.txt, without docking.

`qualification_set.py` writes both files itself after every complex it docks; this script is for reports docked
elsewhere, or for redoing the rates:

    summarize_results.py DIR             every <ligand>_docked_success_rate.csv under DIR (flat, or in <pdb>/ folders:
                                         the reports of parallel CI jobs or of several machines gathered in one place)
                                         -> DIR/qualification_summary.csv and DIR/success_rate.txt.
                                         With --config-dir, the complexes are counted against its TOMLs and the
                                         missing ones listed.
    summarize_results.py --from-summary CSV
                                         only success_rate.txt, in the current directory, from an existing
                                         qualification_summary.csv (e.g. several joined by hand).

success_rate.txt holds, for every RMSD cutoff 0.5, 1.0, ... up to --rmsd-cutoff, the % of complexes whose topK_RMSD
is within it, for every topK, and whose overall_min_RMSD is within it (column overall).

Examples:
    python3 scripts/summarize_results.py ci_collect/steadystate_symmetric --config-dir v2013-coreset/config
    python3 scripts/summarize_results.py gathered_reports --summary-csv results/qualification_summary.csv
    python3 scripts/summarize_results.py --from-summary steadystate/qualification_summary.csv --rmsd-cutoff 3
"""

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from coreset import summary  # noqa: E402


def positive(text: str) -> float:
    value = float(text)
    if not value > 0:
        raise argparse.ArgumentTypeError(f"must be positive, got {text}")
    return value


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Merge DockThor success-rate reports into qualification_summary.csv and success_rate.txt.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "reports", type=Path, nargs="?", metavar="DIR", help="Folder holding <ligand>_docked_success_rate.csv files."
    )
    source.add_argument(
        "--from-summary",
        type=Path,
        metavar="CSV",
        help="Only write success_rate.txt, in the current directory, from this qualification_summary.csv.",
    )
    parser.add_argument(
        "--config-dir",
        type=Path,
        metavar="DIR",
        help="With DIR: count the complexes against the <pdb>.toml of this folder and list the missing ones.",
    )
    parser.add_argument(
        "--summary-csv",
        type=Path,
        metavar="FILE",
        help="With DIR: the summary to write (default DIR/qualification_summary.csv); success_rate.txt goes beside it.",
    )
    parser.add_argument(
        "--rmsd-cutoff",
        type=positive,
        default=3.0,
        metavar="RMSD",
        help="Largest RMSD cutoff of success_rate.txt (default 3.0).",
    )
    args = parser.parse_args(argv)
    if args.from_summary is not None and (args.config_dir is not None or args.summary_csv is not None):
        parser.error("--config-dir and --summary-csv go with DIR, not with --from-summary")
    return args


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.from_summary is not None:
        return summary.generate_summary(args.from_summary, args.rmsd_cutoff)
    return summary.collect(args.reports, args.config_dir, args.summary_csv, args.rmsd_cutoff)


if __name__ == "__main__":
    sys.exit(main())
