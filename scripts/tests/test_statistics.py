"""coreset.statistics: the scoring power of a run -- correlation and error of the score against the affinity."""

import csv
import math

import pytest

from coreset.statistics import (
    MIN_POINTS,
    STATISTIC_COLUMNS,
    least_squares,
    model_of,
    pearson,
    ranks,
    rmse,
    scoring_power_lines,
    scoring_power_table,
    spearman,
    statistics,
    write_scoring_power_csv,
    write_scoring_power_txt,
)

SCORES = ["docktscore_general_dock", "docktscore_ppi"]


def test_pearson_of_a_perfect_line():
    assert pearson([1.0, 2.0, 3.0], [2.0, 4.0, 6.0]) == pytest.approx(1.0)
    assert pearson([1.0, 2.0, 3.0], [-2.0, -4.0, -6.0]) == pytest.approx(-1.0)


def test_pearson_of_a_known_series():
    """A worked example: r of (1,2,3,4) against (2,4,5,4) is 3.5 / sqrt(5 * 4.75) = 0.7181848."""
    assert pearson([1.0, 2.0, 3.0, 4.0], [2.0, 4.0, 5.0, 4.0]) == pytest.approx(0.7181848, abs=1e-6)


def test_a_constant_series_has_no_correlation():
    assert pearson([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]) is None
    assert pearson([1.0, 2.0, 3.0], [1.0, 1.0, 1.0]) is None
    assert least_squares([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]) is None


def test_a_constant_series_is_seen_as_constant_not_as_rounding_noise():
    """The deviations of three equal values from their own mean are ~1e-16, not 0, which used to give a tiny
    nonsense correlation (and so a red badge) instead of `undefined`."""
    same = -6.946457969891734
    result = statistics([same, same, same], [-4.5, -8.0, -8.7])
    assert (result["pearson"], result["r2"], result["spearman"], result["rmse_fit"]) == (None, None, None, None)


def test_too_few_points_leave_the_statistics_undefined():
    """Two points always fit a line perfectly, so a correlation of them would say nothing."""
    assert len([1.0, 2.0]) < MIN_POINTS
    assert pearson([1.0, 2.0], [1.0, 2.0]) is None
    result = statistics([1.0, 2.0], [1.0, 2.0])
    assert (result["n"], result["pearson"], result["rmse_fit"]) == (2, None, None)
    assert result["rmse_raw"] == pytest.approx(0.0)  # the plain error needs no fit


def test_ranks_share_the_average_rank_on_a_tie():
    assert ranks([10.0, 30.0, 20.0]) == [1.0, 3.0, 2.0]
    assert ranks([10.0, 20.0, 20.0, 40.0]) == [1.0, 2.5, 2.5, 4.0]


def test_spearman_sees_a_monotonic_relation_pearson_misses():
    """Ranks only: a curve that always rises correlates perfectly, however bent."""
    xs = [1.0, 2.0, 3.0, 4.0]
    ys = [1.0, 8.0, 27.0, 64.0]
    assert spearman(xs, ys) == pytest.approx(1.0)
    assert pearson(xs, ys) < 1.0


def test_least_squares_and_the_two_errors():
    """The fit takes the score to the measured value, so its RMSE is smaller than reading the score as the value."""
    scores = [-10.0, -8.0, -6.0]
    measured = [-12.0, -10.0, -8.0]  # measured = score - 2
    slope, intercept = least_squares(scores, measured)
    assert (slope, intercept) == (pytest.approx(1.0), pytest.approx(-2.0))
    result = statistics(scores, measured)
    assert result["rmse_fit"] == pytest.approx(0.0)
    assert result["rmse_raw"] == pytest.approx(2.0)
    assert rmse([1.0, 2.0], [1.0, 4.0]) == pytest.approx(math.sqrt(2.0))  # sqrt((0 + 4) / 2)


def rows():
    """Summary rows as the CSV holds them: text, one complex without a score and one without an affinity. The two
    models are on different scales, as general_dock and ppi are."""
    return [
        {"complex": "1bcu", "docktscore_general_dock": "-7.0", "docktscore_ppi": "-14.0", "delta_g_kcal_mol": "-4.5"},
        {"complex": "3ehy", "docktscore_general_dock": "-6.9", "docktscore_ppi": "-13.8", "delta_g_kcal_mol": "-8.0"},
        {"complex": "10gs", "docktscore_general_dock": "-4.4", "docktscore_ppi": "-8.8", "delta_g_kcal_mol": "-8.7"},
        {"complex": "4tmn", "docktscore_general_dock": "", "delta_g_kcal_mol": "-9.0"},  # not rescored
        {"complex": "1ps3", "docktscore_general_dock": "-5.0"},  # no affinity in the index
    ]


def test_model_of_a_score_column():
    assert model_of("docktscore_general_dock") == "general_dock"
    assert model_of("docktscore") == "docktscore"  # a run of before the models were reported separately


def test_scoring_power_counts_only_the_complexes_with_both_values():
    lines = scoring_power_lines(rows(), SCORES, "delta_g_kcal_mol", docked=5)
    assert "Complexes with results: 5" in lines
    assert "Complexes scored: 3" in lines
    assert not any("4tmn" in line or "1ps3" in line for line in lines)


def test_scoring_power_has_a_block_per_model_and_a_table_comparing_them():
    lines = scoring_power_lines(rows(), SCORES, "delta_g_kcal_mol", docked=3)
    text = "\n".join(lines)
    assert "## general_dock (docktscore_general_dock)" in text
    assert "## ppi (docktscore_ppi)" in text
    for label in ("Pearson r", "R^2", "Spearman rho", "RMSE (fit)", "RMSE (raw)", "Fit"):
        assert label in text
    # The comparison table has a line per model, before the blocks.
    table = [line for line in lines if line.strip().startswith(("general_dock", "ppi"))]
    assert len(table) == 2
    assert "N/A" not in text


def test_the_same_scale_free_correlation_for_a_model_on_another_scale():
    """ppi here is general_dock doubled: the correlations are identical, the raw RMSE is not."""
    table = {row["model"]: row for row in scoring_power_table(rows(), SCORES, "delta_g_kcal_mol")}
    assert table["ppi"]["pearson"] == pytest.approx(table["general_dock"]["pearson"])
    assert table["ppi"]["spearman"] == pytest.approx(table["general_dock"]["spearman"])
    assert table["ppi"]["rmse_raw"] > table["general_dock"]["rmse_raw"]
    assert table["ppi"]["rmse_fit"] == pytest.approx(table["general_dock"]["rmse_fit"])


def test_scoring_power_says_so_when_there_is_nothing_to_correlate():
    """One complex: every statistic is N/A and the file says why, rather than printing a number."""
    text = "\n".join(scoring_power_lines(rows()[:1], SCORES, "delta_g_kcal_mol", docked=1))
    assert "Pearson r      : N/A" in text
    assert f"fewer than {MIN_POINTS} complexes" in text


def test_write_scoring_power_txt(tmp_path):
    path = tmp_path / "scoring_power.txt"
    write_scoring_power_txt(path, rows(), SCORES, "delta_g_kcal_mol", 5)
    assert path.read_text().startswith("Scoring power: the DockTScore of the best-scored pose against")
    assert path.read_text().endswith("\n")


def test_write_scoring_power_csv(tmp_path):
    """The machine-readable table the badges are made from: one row per model, an undefined statistic empty."""
    path = tmp_path / "scoring_power.csv"
    write_scoring_power_csv(path, rows(), SCORES, "delta_g_kcal_mol")
    with open(path, newline="") as stream:
        table = list(csv.DictReader(stream))
    assert [row["model"] for row in table] == ["general_dock", "ppi"]
    assert list(table[0]) == STATISTIC_COLUMNS
    assert table[0]["n"] == "3"
    assert float(table[0]["pearson"]) == pytest.approx(-0.6547, abs=1e-3)

    write_scoring_power_csv(path, rows()[:1], SCORES, "delta_g_kcal_mol")
    with open(path, newline="") as stream:
        assert list(csv.DictReader(stream))[0]["pearson"] == ""
