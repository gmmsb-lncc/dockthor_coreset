"""scripts/rescore_results.py and scripts/ci/badge.py: rescoring a finished run, and the badges of its statistics."""

import json
import sys
from pathlib import Path

import pytest

import rescore_results
from coreset.scoring import DEFAULT_MODELS, column_of
from coreset.summary import read_summary
from test_scoring import POSES, TOPOLOGY, fake_docktscore

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "ci"))
import badge  # noqa: E402

AFFINITIES = "complex,log_affinity,delta_g_kcal_mol\n1abc,5.85,-7.979\n2xyz,4.0,-5.456\n3def,6.4,-8.73\n"


def make_run(tmp_path, complexes=("1abc", "2xyz", "3def")):
    """A set of the complexes and a finished run of them: poses and a success-rate report each."""
    dataset = tmp_path / "set"
    docked = tmp_path / "runs"
    for index, pdb in enumerate(complexes):
        (dataset / pdb / "receptors").mkdir(parents=True)
        (dataset / pdb / "receptors" / f"{pdb}.in").write_text(TOPOLOGY)
        out = docked / pdb
        out.mkdir(parents=True)
        (out / f"{pdb}_ligand_docked.mol2").write_text(POSES)
        (out / f"{pdb}_ligand_docked_success_rate.csv").write_text(
            f"complex,execution_time(s),num_top,top1_E,top1_RMSD\n{pdb},1,1,-9,{1.0 + index}\n"
        )
    (dataset / "index").mkdir(parents=True)
    (dataset / "index" / "affinity.csv").write_text(AFFINITIES)
    return dataset, docked


def run(argv):
    return rescore_results.main(argv)


def test_rescore_results_adds_the_scores_and_the_statistics(tmp_path, capsys):
    dataset, docked = make_run(tmp_path)
    docktscore = fake_docktscore(tmp_path / "docktscore.py")

    assert run([str(docked), "--dataset", str(dataset), "--docktscore-path", str(docktscore)]) == 0

    rows = read_summary(docked / "qualification_summary.csv")
    assert [row["complex"] for row in rows] == ["1abc", "2xyz", "3def"]
    for row in rows:
        # Everything DockThor reported is kept, and every model is added.
        assert row["top1_RMSD"]
        assert [row[column_of(model)] for model in DEFAULT_MODELS] == ["-6.946457969891734"] * len(DEFAULT_MODELS)
        assert row["delta_g_kcal_mol"]
    assert (docked / "success_rate.txt").is_file()
    power = (docked / "scoring_power.txt").read_text()
    for model in DEFAULT_MODELS:
        assert f"## {model} ({column_of(model)})" in power
    table = read_summary(docked / "scoring_power.csv")
    assert [row["model"] for row in table] == DEFAULT_MODELS
    assert table[0]["n"] == "3"
    # A constant score cannot correlate with anything, so the statistics are empty rather than wrong.
    assert table[0]["pearson"] == ""


def test_rescore_results_takes_the_models_and_the_complexes_given(tmp_path, capsys):
    dataset, docked = make_run(tmp_path)
    docktscore = fake_docktscore(tmp_path / "docktscore.py")
    args = [
        str(docked),
        "--dataset",
        str(dataset),
        "--docktscore-path",
        str(docktscore),
        "--docktscore-models",
        "protease",
        "--complex",
        "1abc",
    ]

    assert run(args) == 0
    rows = {row["complex"]: row for row in read_summary(docked / "qualification_summary.csv")}
    assert rows["1abc"]["docktscore_protease"] == "-6.946457969891734"
    assert rows["2xyz"]["docktscore_protease"] == ""  # not asked for, so not scored
    assert "docktscore_general_dock" not in rows["1abc"]
    assert [row["model"] for row in read_summary(docked / "scoring_power.csv")] == ["protease"]


def test_rescore_results_skips_what_is_already_scored(tmp_path, capsys):
    dataset, docked = make_run(tmp_path)
    docktscore = fake_docktscore(tmp_path / "docktscore.py")
    args = [str(docked), "--dataset", str(dataset), "--docktscore-path", str(docktscore), "--skip-done"]
    assert run(args) == 0
    capsys.readouterr()

    assert run(args) == 0  # the scores are on disk in <pdb>_scores.csv, so the second run has nothing to do
    out = capsys.readouterr().out
    assert out.count("skip, already scored") == 3
    assert "Rescored 0 of 3" in out


def test_rescore_results_needs_a_run_with_poses_and_a_set(tmp_path, capsys):
    dataset, docked = make_run(tmp_path)
    docktscore = fake_docktscore(tmp_path / "docktscore.py")
    with pytest.raises(SystemExit):  # neither --dataset nor --config-dir
        run([str(docked), "--docktscore-path", str(docktscore)])
    with pytest.raises(SystemExit):  # the binary must exist
        run([str(docked), "--dataset", str(dataset), "--docktscore-path", str(tmp_path / "absent")])

    empty = tmp_path / "empty"
    empty.mkdir()
    assert run([str(empty), "--dataset", str(dataset), "--docktscore-path", str(docktscore)]) == 1
    assert "no *_docked_success_rate.csv" in capsys.readouterr().err


def test_badges_of_the_scoring_power(tmp_path, capsys):
    """One badge per model and statistic, with the color of the value and the version as the label."""
    summary = tmp_path / "qualification_summary.csv"
    summary.write_text(
        "complex,overall_min_RMSD,top1_RMSD\n1abc,0.5,0.5\n2xyz,1.0,3.0\n"  # 100% overall, 50% top1 at 2 Å
    )
    power = tmp_path / "scoring_power.csv"
    power.write_text(
        "model,n,pearson,r2,spearman,rmse_fit,rmse_raw,slope,intercept\n"
        "general_dock,195,0.72,0.5184,0.70,1.44,2.0,1.0,0.0\n"
        "ppi,195,-0.35,0.1225,-0.30,3.8,9.0,1.0,0.0\n"
        "protease,2,,,,,,,\n"  # too few complexes: no statistic, so no badge
    )
    prefix = tmp_path / "badges" / "steadystate_symmetric"

    assert badge.main([str(summary), "--version", "2.6.0", "-o", str(prefix), "--scoring-power", str(power)]) == 0

    def read(name):
        return json.loads(prefix.with_name(f"{prefix.name}_{name}.json").read_text())

    assert read("top1")["message"] == "50.0% ≤ 2 Å"
    assert read("general_dock_pearson") == {
        "schemaVersion": 1,
        "label": "2.6.0",
        "message": "r = 0.72",
        "color": "brightgreen",
    }
    assert read("general_dock_r2")["message"] == "R² = 0.52"
    assert read("general_dock_spearman")["message"] == "ρ = 0.70"
    assert read("general_dock_rmse") == {
        "schemaVersion": 1,
        "label": "2.6.0",
        "message": "1.44 kcal/mol",
        "color": "brightgreen",
    }
    # A correlation in the wrong direction is red, and a large error too.
    assert read("ppi_pearson")["color"] == "red"
    assert read("ppi_rmse")["color"] == "red"
    assert not prefix.with_name(f"{prefix.name}_protease_pearson.json").exists()


def test_badges_need_a_scoring_power_with_a_statistic(tmp_path, capsys):
    summary = tmp_path / "qualification_summary.csv"
    summary.write_text("complex,overall_min_RMSD,top1_RMSD\n1abc,0.5,0.5\n")
    power = tmp_path / "scoring_power.csv"
    power.write_text("model,n,pearson,r2,spearman,rmse_fit,rmse_raw,slope,intercept\nppi,2,,,,,,,\n")
    prefix = tmp_path / "badges" / "setup"

    assert badge.main([str(summary), "--version", "2.6.0", "-o", str(prefix), "--scoring-power", str(power)]) == 1
    assert "no statistic to make a badge from" in capsys.readouterr().err
    assert badge.main([str(summary), "--version", "2.6.0", "-o", str(prefix), "--scoring-power", str(tmp_path / "x")])
