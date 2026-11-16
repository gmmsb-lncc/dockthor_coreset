#!/usr/bin/env python3
"""Write the shields.io endpoint badges of one setup: the success rates, and the scoring power of each model.

Each badge is one JSON (an endpoint JSON holds one badge). The **success rates** come from
qualification_summary.csv, at one RMSD cutoff (default 2.0 Å), as in success_rate.txt:

- <prefix>_overall.json: the % of complexes whose overall_min_RMSD (the lowest-RMSD leader of all the leaders, so of
  any of the runs) is within the cutoff, the `overall` column;
- <prefix>_top1.json: the % of complexes whose top1_RMSD (the best-scored pose) is within the cutoff, the `top1`
  column: the usual docking success rate.

With `--scoring-power scoring_power.csv` (written by the rescoring, README 5.5), the **scoring power** of every
DockTScore model is added, four badges per model -- `<prefix>_<model>_pearson.json`, `_r2.json`, `_spearman.json`
and `_rmse.json` -- which is one row per model and one column per statistic in a README table.

The label is the DockThor version, the message the rate or the statistic. Displayed with:

    https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/gmmsb-lncc/dockthor_coreset/main/dockthor_results/badges/<setup>_overall.json

Example:
    python3 scripts/ci/badge.py dockthor_results/2.2.1/steadystate_symmetric_qualification_summary.csv \\
        --version 2.2.1 -o dockthor_results/badges/steadystate_symmetric \\
        --scoring-power dockthor_results/2.2.1/steadystate_symmetric_scoring_power.csv
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from coreset.summary import read_summary  # noqa: E402

# (lowest rate in %, shields.io color), first match wins
COLORS = [(90.0, "brightgreen"), (80.0, "green"), (70.0, "yellowgreen"), (60.0, "yellow"), (50.0, "orange")]
# A correlation: the signed value, so an anti-correlated model (the wrong direction) is red.
CORRELATION_COLORS = [(0.7, "brightgreen"), (0.6, "green"), (0.5, "yellowgreen"), (0.4, "yellow"), (0.3, "orange")]
# An error in kcal/mol: lower is better, so the thresholds are ceilings.
RMSE_COLORS = [(1.5, "brightgreen"), (2.0, "green"), (2.5, "yellowgreen"), (3.0, "yellow"), (3.5, "orange")]

# badge name: column of qualification_summary.csv
METRICS = {"overall": "overall_min_RMSD", "top1": "top1_RMSD"}
# badge name: (column of scoring_power.csv, how the message reads)
STATISTICS = {
    "pearson": ("pearson", "r = {:.2f}"),
    "r2": ("r2", "R² = {:.2f}"),
    "spearman": ("spearman", "ρ = {:.2f}"),
    "rmse": ("rmse_fit", "{:.2f} kcal/mol"),
}


def color_of(rate: float) -> str:
    return next((color for floor, color in COLORS if rate >= floor), "red")


def correlation_color(value: float) -> str:
    return next((color for floor, color in CORRELATION_COLORS if value >= floor), "red")


def rmse_color(value: float) -> str:
    return next((color for ceiling, color in RMSE_COLORS if value <= ceiling), "red")


def write_badge(path: Path, label: str, message: str, color: str) -> None:
    badge = {"schemaVersion": 1, "label": label, "message": message, "color": color}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(badge, ensure_ascii=False, indent=2) + "\n")


def scoring_power_badges(scoring_power: Path, version: str, prefix: Path) -> int:
    """One badge per model and statistic of a scoring_power.csv; the number written."""
    written = 0
    for row in read_summary(scoring_power):
        model = row.get("model", "")
        for name, (column, template) in STATISTICS.items():
            if not row.get(column):  # a statistic that is not defined (fewer than three complexes, say)
                continue
            value = float(row[column])
            color = rmse_color(value) if name == "rmse" else correlation_color(value)
            output = prefix.with_name(f"{prefix.name}_{model}_{name}.json")
            write_badge(output, version, template.format(value), color)
            print(f"{output}: {model} {name} {version}: {template.format(value)}")
            written += 1
    return written


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("summary", type=Path, help="qualification_summary.csv")
    parser.add_argument("--version", required=True, help="DockThor version, the badge's label")
    parser.add_argument("--cutoff", type=float, default=2.0, help="RMSD cutoff in Å")
    parser.add_argument(
        "-o", "--output-prefix", type=Path, required=True, help="Badges to write: <prefix>_overall.json, _top1.json"
    )
    parser.add_argument(
        "--scoring-power",
        type=Path,
        metavar="CSV",
        help="scoring_power.csv: also write <prefix>_<model>_{pearson,r2,spearman,rmse}.json",
    )
    args = parser.parse_args(argv)

    reports = read_summary(args.summary)
    for name, column in METRICS.items():
        values = [float(r[column]) for r in reports if r.get(column)]
        if not values:
            print(f"[error] {args.summary}: no {column} values", file=sys.stderr)
            return 1
        rate = 100 * sum(v <= args.cutoff for v in values) / len(values)
        output = args.output_prefix.with_name(f"{args.output_prefix.name}_{name}.json")
        write_badge(output, args.version, f"{rate:.1f}% ≤ {args.cutoff:g} Å", color_of(rate))
        print(f"{output}: {name} {args.version}: {rate:.1f}% ≤ {args.cutoff:g} Å ({len(values)} complexes)")
    if args.scoring_power is not None:
        if not args.scoring_power.is_file():
            print(f"[error] {args.scoring_power}: not found", file=sys.stderr)
            return 1
        if not scoring_power_badges(args.scoring_power, args.version, args.output_prefix):
            print(f"[error] {args.scoring_power}: no statistic to make a badge from", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
