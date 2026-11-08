#!/usr/bin/env python3
"""Write the shields.io endpoint badges of one setup from its qualification_summary.csv.

Two badges, one JSON each (an endpoint JSON holds one badge), with the success rate at one RMSD cutoff (default
2.0 Å), as in success_rate.txt:

- <prefix>_overall.json: the % of complexes whose overall_min_RMSD (the lowest-RMSD leader of all the leaders, so of
  any of the runs) is within the cutoff, the `overall` column;
- <prefix>_top1.json: the % of complexes whose top1_RMSD (the best-scored pose) is within the cutoff, the `top1`
  column: the usual docking success rate.

The label is the DockThor version, the message the rate. Displayed with:

    https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/gmmsb-lncc/dockthor_coreset/main/dockthor_results/badges/<setup>_overall.json

Example:
    python3 scripts/ci/badge.py dockthor_results/2.2.1/steadystate_symmetric_qualification_summary.csv \\
        --version 2.2.1 -o dockthor_results/badges/steadystate_symmetric
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import qualification_set as qs  # noqa: E402  (scripts/qualification_set.py)

# (lowest rate in %, shields.io color), first match wins
COLORS = [(90.0, "brightgreen"), (80.0, "green"), (70.0, "yellowgreen"), (60.0, "yellow"), (50.0, "orange")]

# badge name: column of qualification_summary.csv
METRICS = {"overall": "overall_min_RMSD", "top1": "top1_RMSD"}


def color_of(rate: float) -> str:
    return next((color for floor, color in COLORS if rate >= floor), "red")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("summary", type=Path, help="qualification_summary.csv")
    parser.add_argument("--version", required=True, help="DockThor version, the badge's label")
    parser.add_argument("--cutoff", type=float, default=2.0, help="RMSD cutoff in Å")
    parser.add_argument(
        "-o", "--output-prefix", type=Path, required=True, help="Badges to write: <prefix>_overall.json, _top1.json"
    )
    args = parser.parse_args(argv)

    reports = qs.read_summary(args.summary)
    for name, column in METRICS.items():
        values = [float(r[column]) for r in reports if r.get(column)]
        if not values:
            print(f"[error] {args.summary}: no {column} values", file=sys.stderr)
            return 1
        rate = 100 * sum(v <= args.cutoff for v in values) / len(values)
        badge = {
            "schemaVersion": 1,
            "label": args.version,
            "message": f"{rate:.1f}% ≤ {args.cutoff:g} Å",
            "color": color_of(rate),
        }
        output = args.output_prefix.with_name(f"{args.output_prefix.name}_{name}.json")
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(badge, ensure_ascii=False, indent=2) + "\n")
        print(f"{output}: {name} {badge['label']}: {badge['message']} ({len(values)} complexes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
