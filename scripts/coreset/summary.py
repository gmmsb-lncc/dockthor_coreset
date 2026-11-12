"""The success-rate reports: DockThor's <ligand>_docked_success_rate.csv, their consolidation qualification_summary.csv,
and success_rate.txt; the --generate-summary and --collect commands."""

import csv
import re
import sys
from pathlib import Path

SUCCESS_RATE_CSV = (
    "_docked_success_rate.csv"  # suffix of DockThor's --success-rate report, <ligand>_docked_success_rate.csv
)
OVERALL = "overall"  # success-rate column of overall_min_RMSD, before the topK ones


def read_success_rate_file(path: Path) -> dict[str, str]:
    """{column: value} of one <ligand>_docked_success_rate.csv; empty if it is malformed."""
    with open(path, newline="") as stream:
        rows = [row for row in csv.reader(stream) if row]
    return dict(zip(rows[0], rows[1], strict=False))  # a short row: the columns it has if len(rows) >= 2 else {}


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


def read_summary(path: Path) -> list[dict]:
    """The rows of a qualification_summary.csv, one {column: value} per complex."""
    with open(path, newline="") as stream:
        return [row for row in csv.DictReader(stream) if any(row.values())]


def success_rate_cutoffs(max_rmsd: float) -> list[float]:
    """0.5, 1.0, ... up to max_rmsd rounded up to the 0.5 grid (0.3 -> [0.5], 2.2 -> [..., 2.5]),
    the cutoffs of DockThor's --success-rate-by-run."""
    cutoffs = []
    while not cutoffs or cutoffs[-1] < max_rmsd - 1e-4:
        cutoffs.append(0.5 * (len(cutoffs) + 1))
    return cutoffs


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


def generate_summary(summary_csv: Path, rmsd_cutoff: float) -> int:
    """--generate-summary: success_rate.txt, in the current directory, from an existing summary CSV."""
    summary_path = summary_csv.expanduser().resolve()
    try:
        reports = read_summary(summary_path)
    except OSError as error:
        print(f"[error] {error}", file=sys.stderr)
        return 1
    if not top_ks(reports):
        print(f"[error] {summary_path}: no topK_RMSD column, not a qualification_summary.csv", file=sys.stderr)
        return 1
    success_rate_path = Path.cwd() / "success_rate.txt"
    write_success_rate_txt(success_rate_path, reports, success_rate_cutoffs(rmsd_cutoff), len(reports))
    print(f"Success rates of {len(reports)} complex(es) from {summary_path}: {success_rate_path}")
    return 0


def collect_reports(root: Path) -> tuple[list[dict], set[str]]:
    """One report per complex from every <ligand>_docked_success_rate.csv under root, sorted by complex, and the
    names of those complexes (the report's complex column, else the file name). A malformed file or a second
    report of a complex is skipped with a warning."""
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
    return reports, seen


def collect(directory: Path, config_dir: Path | None, summary_csv: Path | None, rmsd_cutoff: float) -> int:
    """--collect DIR: qualification_summary.csv and success_rate.txt from every <ligand>_docked_success_rate.csv
    under DIR (flat, or in <pdb>/ folders), e.g. the per-complex reports of parallel runs gathered into one folder.
    With config_dir, the complexes are counted against its TOMLs and the missing ones are listed."""
    root = directory.expanduser().resolve()
    if not root.is_dir():
        print(f"[error] --collect {root} is not a directory", file=sys.stderr)
        return 1

    reports, seen = collect_reports(root)
    if not reports or not top_ks(reports):
        print(f"[error] no *{SUCCESS_RATE_CSV} with topK_RMSD columns under {root}", file=sys.stderr)
        return 1

    expected = []
    if config_dir is not None:
        toml_dir = config_dir.expanduser().resolve()
        if not toml_dir.is_dir():
            print(f"[error] --config-dir {toml_dir} is not a directory", file=sys.stderr)
            return 1
        expected = sorted(p.stem for p in toml_dir.glob("*.toml"))
    num_run = len(expected) if expected else len(reports)

    summary_path = summary_csv.expanduser().resolve() if summary_csv else root / "qualification_summary.csv"
    success_rate_path = summary_path.with_name("success_rate.txt")
    write_summary(summary_path, reports)
    write_success_rate_txt(success_rate_path, reports, success_rate_cutoffs(rmsd_cutoff), num_run)
    print(f"Collected {len(reports)} of {num_run} complex(es) from {root}")
    print(f"Summary: {summary_path}")
    print(f"Success rates: {success_rate_path}")
    if missing := [pdb for pdb in expected if pdb not in seen]:
        print(f"[warn] {len(missing)} complex(es) without results: {' '.join(missing)}")
    return 0
