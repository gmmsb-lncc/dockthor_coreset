"""The scoring power of a qualification run: how well each DockTScore model reproduces the measured affinities.

The docking power of a set is its success rate, the % of complexes whose pose is within an RMSD cutoff
(`success_rate.txt`, `coreset/summary.py`). Its **scoring power** is the agreement between the score of the
best-scored pose and the measured binding affinity, over the complexes that have both, and is what
`scoring_power.txt` and `scoring_power.csv` report, **one row per DockTScore model**:

| Statistic | What it says |
| --------- | ------------ |
| Pearson r | the linear correlation of the score with the measured value; CASF's scoring power |
| R^2 | r^2, the share of the variance of the measured values the score explains |
| Spearman rho | the same correlation on ranks: insensitive to how the two scales differ |
| RMSE (fit) | the error in kcal/mol after the least-squares line score -> measured value |
| RMSE (raw) | the error of reading the score as the affinity directly, with no fit |

The score is compared with `delta_g_kcal_mol`, the measured affinity in the same unit, so the raw RMSE means
something; `|r|` is the same against `log_affinity`, which is only that value divided by a constant. A model fitted
on another scale (`protease`, `ppi`) is told apart by its raw RMSE, not by its correlation, which is scale-free.

Everything here is the plain stdlib, and n < 3 leaves the statistics undefined rather than reporting a number from
two points.
"""

import csv
import math

#: Fewest complexes a correlation is computed from; below it the statistics are reported as N/A.
MIN_POINTS = 3
#: Columns of scoring_power.csv, and the statistics of the comparison table, in order.
STATISTIC_COLUMNS = ["model", "n", "pearson", "r2", "spearman", "rmse_fit", "rmse_raw", "slope", "intercept"]


def is_constant(values: list[float]) -> bool:
    """True when every value is the same: such a series has no correlation and no least-squares line. Checked on
    the values themselves, because the deviations from their mean are rounding noise rather than zero."""
    return min(values) == max(values)


def pearson(xs: list[float], ys: list[float]) -> float | None:
    """Pearson correlation of two equally long series; None when either is constant or there are too few points."""
    if len(xs) < MIN_POINTS or is_constant(xs) or is_constant(ys):
        return None
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    dx = [x - mean_x for x in xs]
    dy = [y - mean_y for y in ys]
    spread = math.sqrt(sum(value * value for value in dx) * sum(value * value for value in dy))
    return sum(a * b for a, b in zip(dx, dy, strict=True)) / spread if spread > 0 else None


def ranks(values: list[float]) -> list[float]:
    """The ranks of the values, ties sharing their average rank (Spearman's ranking)."""
    order = sorted(range(len(values)), key=lambda index: values[index])
    result = [0.0] * len(values)
    start = 0
    while start < len(order):
        stop = start
        while stop + 1 < len(order) and values[order[stop + 1]] == values[order[start]]:
            stop += 1
        average = (start + stop) / 2 + 1
        for index in order[start : stop + 1]:
            result[index] = average
        start = stop + 1
    return result


def spearman(xs: list[float], ys: list[float]) -> float | None:
    """Spearman rank correlation, the Pearson correlation of the ranks."""
    return pearson(ranks(xs), ranks(ys)) if len(xs) >= MIN_POINTS else None


def least_squares(xs: list[float], ys: list[float]) -> tuple[float, float] | None:
    """(slope, intercept) of the least-squares line y = slope * x + intercept; None when x is constant."""
    if len(xs) < MIN_POINTS or is_constant(xs):
        return None
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    variance = sum((x - mean_x) ** 2 for x in xs)
    if variance <= 0:
        return None
    slope = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys, strict=True)) / variance
    return slope, mean_y - slope * mean_x


def rmse(predicted: list[float], measured: list[float]) -> float:
    """Root-mean-square error of two equally long series."""
    squares = [(p - m) ** 2 for p, m in zip(predicted, measured, strict=True)]
    return math.sqrt(sum(squares) / len(squares))


def statistics(scores: list[float], measured: list[float]) -> dict:
    """{n, pearson, r2, spearman, slope, intercept, rmse_fit, rmse_raw} of a score series and the measured values.

    rmse_fit is the error after the least-squares line from the score to the measured value, rmse_raw the error of
    using the score as the measured value directly. Correlations are None when they are undefined."""
    result = {"n": len(scores), "pearson": None, "r2": None, "spearman": None}
    result |= {"slope": None, "intercept": None, "rmse_fit": None, "rmse_raw": None}
    if not scores:
        return result
    result["rmse_raw"] = rmse(scores, measured)
    if (r := pearson(scores, measured)) is not None:
        result["pearson"] = r
        result["r2"] = r * r
    result["spearman"] = spearman(scores, measured)
    if (line := least_squares(scores, measured)) is not None:
        slope, intercept = line
        result |= {"slope": slope, "intercept": intercept}
        result["rmse_fit"] = rmse([slope * score + intercept for score in scores], measured)
    return result


def points_of(rows: list[dict], score_column: str, measured_column: str) -> list[tuple[str, float, float]]:
    """(complex, score, measured value) of every row holding both numbers, in the order of the rows."""
    points = []
    for row in rows:
        try:
            points.append((row.get("complex", ""), float(row[score_column]), float(row[measured_column])))
        except (KeyError, TypeError, ValueError):
            continue
    return points


def model_statistics(rows: list[dict], score_column: str, measured_column: str) -> dict:
    """The statistics of one score column against the measured one, with its points under `points`."""
    points = points_of(rows, score_column, measured_column)
    result = statistics([score for _, score, _ in points], [value for _, _, value in points])
    return result | {"points": points}


def number(value: float | None, digits: int = 4) -> str:
    """A statistic for the report: N/A when it is undefined."""
    return "N/A" if value is None else f"{value:.{digits}f}"


def model_of(score_column: str) -> str:
    """The DockTScore model a score column belongs to: `docktscore_general_dock` -> `general_dock`."""
    return score_column.split("_", 1)[1] if "_" in score_column else score_column


def scoring_power_table(rows: list[dict], score_columns: list[str], measured_column: str) -> list[dict]:
    """One {model, n, pearson, ...} per score column: the rows of scoring_power.csv, in the order given.

    The values are floats or None, so a reader does not parse N/A; `write_scoring_power_csv` prints them."""
    table = []
    for column in score_columns:
        result = model_statistics(rows, column, measured_column)
        table.append({"model": model_of(column)} | {key: result[key] for key in STATISTIC_COLUMNS[1:]})
    return table


def write_scoring_power_csv(path, rows: list[dict], score_columns: list[str], measured_column: str) -> None:
    """scoring_power.csv: one row per model, the columns of STATISTIC_COLUMNS, an undefined statistic left empty."""
    with open(path, "w", newline="") as stream:
        writer = csv.DictWriter(stream, STATISTIC_COLUMNS)
        writer.writeheader()
        for row in scoring_power_table(rows, score_columns, measured_column):
            writer.writerow({key: "" if row[key] is None else row[key] for key in STATISTIC_COLUMNS})


def comparison_lines(table: list[dict]) -> list[str]:
    """The table that compares the models: one line each, the statistics as columns."""
    header = f"{'model':>14}{'n':>6}{'Pearson r':>12}{'R^2':>9}{'Spearman':>10}{'RMSE fit':>10}{'RMSE raw':>10}"
    lines = [header]
    for row in table:
        lines.append(
            f"{row['model']:>14}{row['n']:>6}{number(row['pearson']):>12}{number(row['r2'], 3):>9}"
            f"{number(row['spearman'], 3):>10}{number(row['rmse_fit'], 3):>10}{number(row['rmse_raw'], 3):>10}"
        )
    return lines


def model_lines(rows: list[dict], score_column: str, measured_column: str) -> list[str]:
    """One model's block: its statistics, its fit, and one line per complex with the fitted value and the residual."""
    result = model_statistics(rows, score_column, measured_column)
    points = result["points"]
    lines = [
        f"## {model_of(score_column)} ({score_column})",
        f"Complexes scored: {result['n']}",
        f"Pearson r      : {number(result['pearson'])}",
        f"R^2            : {number(result['r2'])}",
        f"Spearman rho   : {number(result['spearman'])}",
        f"RMSE (fit)     : {number(result['rmse_fit'])} kcal/mol",
        f"RMSE (raw)     : {number(result['rmse_raw'])} kcal/mol",
    ]
    if result["slope"] is not None:
        sign = "+" if result["intercept"] >= 0 else "-"
        lines.append(
            f"Fit            : {measured_column} = {result['slope']:.4f} * {score_column} "
            f"{sign} {abs(result['intercept']):.4f}"
        )
    if result["n"] < MIN_POINTS:
        lines.append(f"[warn] fewer than {MIN_POINTS} complexes: the correlations are not defined")
    if not points:
        return lines
    fitted = result["slope"] is not None
    lines += ["", f"{'complex':>8}{'score':>14}{measured_column:>18}{'fitted':>10}{'residual':>10}"]
    for name, score, value in points:
        cells = f"{name:>8}{score:>14.4f}{value:>18.4f}"
        if fitted:
            prediction = result["slope"] * score + result["intercept"]
            cells += f"{prediction:>10.4f}{prediction - value:>10.4f}"
        else:
            cells += f"{'N/A':>10}{'N/A':>10}"
        lines.append(cells)
    return lines


def scoring_power_lines(rows: list[dict], score_columns: list[str], measured_column: str, docked: int) -> list[str]:
    """The lines of scoring_power.txt: the model comparison, then one block per model.

    `rows` are summary rows ({column: text}, `coreset/summary.py`), `docked` the number of complexes with results."""
    lines = [
        f"Scoring power: the DockTScore of the best-scored pose against {measured_column}",
        f"Complexes with results: {docked}",
        "",
        *comparison_lines(scoring_power_table(rows, score_columns, measured_column)),
    ]
    for column in score_columns:
        lines += ["", ""] + model_lines(rows, column, measured_column)
    return lines


def write_scoring_power_txt(path, rows: list[dict], score_columns: list[str], measured_column: str, docked: int):
    """Write scoring_power.txt."""
    with open(path, "w") as stream:
        stream.write("\n".join(scoring_power_lines(rows, score_columns, measured_column, docked)) + "\n")
