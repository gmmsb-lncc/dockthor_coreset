#!/usr/bin/env python3
"""Write a shields.io endpoint badge from a qualification_summary.csv.

The badge shows the overall success rate at one RMSD cutoff (default 2.0 Å): the % of complexes whose
overall_min_RMSD (the lowest-RMSD leader of all the leaders) is within the cutoff, as in the `overall` column of
success_rate.txt. Displayed with:

    https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/gmmsb-lncc/dockthor_coreset/main/dockthor_results/badges/<setup>.json

Example:
    python3 scripts/ci/badge.py dockthor_results/2.1.0/steadystate_symmetric_qualification_summary.csv \\
        --label "SteadyState sym" --version 2.1.0 -o dockthor_results/badges/steadystate_symmetric.json
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import qualification_set as qs  # noqa: E402  (scripts/qualification_set.py)

# (lowest rate in %, shields.io color), first match wins
COLORS = [(90.0, "brightgreen"), (80.0, "green"), (70.0, "yellowgreen"), (60.0, "yellow"), (50.0, "orange")]


def color_of(rate: float) -> str:
    return next((color for floor, color in COLORS if rate >= floor), "red")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("summary", type=Path, help="qualification_summary.csv")
    parser.add_argument("--label", required=True, help="Left side of the badge, e.g. 'SteadyState sym'")
    parser.add_argument("--version", required=True, help="DockThor version, shown on the badge")
    parser.add_argument("--cutoff", type=float, default=2.0, help="RMSD cutoff in Å")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Badge JSON to write")
    args = parser.parse_args(argv)

    reports = qs.read_summary(args.summary)
    values = [float(r["overall_min_RMSD"]) for r in reports if r.get("overall_min_RMSD")]
    if not values:
        print(f"[error] {args.summary}: no overall_min_RMSD values", file=sys.stderr)
        return 1
    rate = 100 * sum(v <= args.cutoff for v in values) / len(values)

    badge = {
        "schemaVersion": 1,
        "label": f"{args.label} {args.version}",
        "message": f"{rate:.1f}% ≤ {args.cutoff:g} Å",
        "color": color_of(rate),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(badge, ensure_ascii=False, indent=2) + "\n")
    print(f"{args.output}: {badge['label']}: {badge['message']} ({len(values)} complexes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
