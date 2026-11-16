"""The layout of a set (one directory per complex) and the three command lines end to end on a small set:
grid_toml_generator.py, qualification_set.py (a dry run) and summarize_results.py."""

import tomllib
from pathlib import Path

import pytest

import grid_toml_generator
import qualification_set
import summarize_results
from coreset.dataset import cofactors_dir, complex_dirs, infer_dataset
from coreset.scoring import DEFAULT_MODELS, column_of
from coreset.summary import read_summary
from test_molecules import CIF, TOP


def make_set(root):
    """Two complexes, one directory each with a folder per role: 1abc with a .top, a cofactor and an extra MOL2, 2xyz
    with only MMFFLigand's .cif and no receptor; plus folders that are not complexes (index/, a run's output)."""
    for directory in ("1abc/receptors", "1abc/ligands", "1abc/cofactors", "2xyz/ligands", "index", "runs/1abc"):
        (root / directory).mkdir(parents=True)
    (root / "1abc/ligands/1abc_ligand.top").write_text(TOP)
    (root / "1abc/ligands/1abc_ligand_mmff.mol2").write_text("@<TRIPOS>MOLECULE\n")
    (root / "1abc/receptors/1abc.in").write_text("receptor\n")
    (root / "1abc/cofactors/1abc_cofactor-01.top").write_text(TOP)
    (root / "1abc/cofactors/1abc_cofactor-01.pdb").write_text("HETATM\n")
    (root / "2xyz/2xyz_ligand_mmff.cif").write_text(CIF)  # outside ligands/: not the ligand of 2xyz
    (root / "2xyz/ligands/2xyz_ligand_mmff.cif").write_text(CIF)
    (root / "runs/1abc/1abc_ligand_docked.csv").write_text("#pose\n")
    return root


def test_complex_dirs(tmp_path):
    root = make_set(tmp_path)
    assert complex_dirs(root) == ["1abc", "2xyz"]
    assert cofactors_dir(root, "1abc") == Path("1abc/cofactors")
    assert cofactors_dir(root, "2xyz") is None


def test_grid_toml_generator(tmp_path, capsys):
    root = make_set(tmp_path)
    assert grid_toml_generator.main([str(root), "--padding", "20"]) == 1  # 2xyz has no receptor
    out, err = capsys.readouterr()
    assert "Skipped unrecognized 1abc/ligands/1abc_ligand_mmff.mol2" in err
    assert "2xyz: missing" in err
    assert not (root / "config/2xyz.toml").exists()

    text = (root / "config/1abc.toml").read_text()
    config = tomllib.loads(text)
    # extent 1.088 x 1.798 x 1.453 + 20: a cube of 21.798 Å, 88 points per axis
    assert config == {
        "output-dir": "vcoreset_run",
        "ligand": str(root.resolve() / "1abc/ligands/1abc_ligand.top"),
        "receptor": str(root.resolve() / "1abc/receptors/1abc.in"),
        "cofactors-path": str(root.resolve() / "1abc/cofactors"),
        "grid-size": [21.798, 21.798, 21.798],
        "grid-center": [14.544, 9.899, 24.2735],
    }  # no grid-rstep: --rstep not given, DockThor's default 0.25 counts the points
    assert text.startswith("## DockThor qualification set: 1abc\n\n")
    assert "681472" in out  # 88^3
    # The ligand and receptor are absolute: the set root of config/ is still the set, not config/ itself.
    assert infer_dataset(root / "config") == root.resolve()

    # An existing TOML is kept; --box-shape ligand and a large padding: over the limit, grid-max-points written.
    assert grid_toml_generator.main([str(root), "--complex", "1abc", "--padding", "40", "--box-shape", "ligand"]) == 0
    assert "exists (use --overwrite)" in capsys.readouterr().out
    assert (
        grid_toml_generator.main(
            [str(root), "--complex", "1abc", "--padding", "40", "--box-shape", "ligand", "--overwrite"]
        )
        == 0
    )
    text = (root / "config/1abc.toml").read_text()
    config = tomllib.loads(text)
    assert config["grid-size"] == [41.088, 41.798, 41.453]
    assert "grid-rstep" not in config
    assert config["grid-max-points"] == 165 * 168 * 166
    assert f"# grid points: {165 * 168 * 166} (165 x 168 x 166), over DockThor's default limit 1000000\n" in text

    # --rstep given: written, and the points counted with it.
    assert (
        grid_toml_generator.main([str(root), "--complex", "1abc", "--padding", "20", "--rstep", "0.25", "--overwrite"])
        == 0
    )
    config = tomllib.loads((root / "config/1abc.toml").read_text())
    assert config["grid-rstep"] == 0.25
    assert "grid-max-points" not in config
    assert (
        grid_toml_generator.main([str(root), "--complex", "1abc", "--padding", "20", "--rstep", "0.2", "--overwrite"])
        == 0
    )
    config = tomllib.loads((root / "config/1abc.toml").read_text())
    assert config["grid-rstep"] == 0.2
    assert config["grid-max-points"] == 109**3  # 21.798 / 0.2 = 108.99 -> 109 points per axis


def test_qualification_set_docks_the_generated_tomls(tmp_path, capsys):
    root = make_set(tmp_path)
    grid_toml_generator.main([str(root)])  # 1abc only: 2xyz has no receptor
    dockthor = tmp_path / "dockthor"
    dockthor.write_text("")
    capsys.readouterr()

    # The set: its config/ by default.
    assert qualification_set.main([str(root), "-o", "runs", "--dockthor-path", str(dockthor), "--dry-run"]) == 0
    out = capsys.readouterr().out
    assert f"[1/1] 1abc: {dockthor} --config {root.resolve()}/config/1abc.toml --output-dir runs --success-rate" in out

    # --config-dir alone: the set root is inferred from the TOMLs' absolute paths, not taken as config/ itself.
    args = ["--config-dir", str(root / "config"), "--dockthor-path", str(dockthor), "--dry-run"]
    assert qualification_set.main(args) == 0
    assert f"Dataset root (from the TOMLs): {root.resolve()}\n" in capsys.readouterr().out

    # No TOMLs yet: a clear error, not an empty run.
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(SystemExit):
        qualification_set.main([str(empty), "--dockthor-path", str(dockthor)])
    assert "grid_toml_generator.py" in capsys.readouterr().err


FAKE_DOCKTHOR = """#!/usr/bin/env python3
# Writes what DockThor writes of one complex: the poses and the --success-rate report.
import sys
from pathlib import Path
config = Path(sys.argv[sys.argv.index("--config") + 1])
out = Path(sys.argv[sys.argv.index("--output-dir") + 1]) / config.stem
out.mkdir(parents=True, exist_ok=True)
(out / f"{config.stem}_ligand_docked.mol2").write_text("@<TRIPOS>MOLECULE\\nbest\\n@<TRIPOS>MOLECULE\\nsecond\\n")
(out / f"{config.stem}_ligand_docked_success_rate.csv").write_text(
    f"complex,execution_time(s),num_top,top1_E,top1_RMSD\\n{config.stem},1,1,-9,1.2\\n"
)
"""

FAKE_DOCKTSCORE = """#!/usr/bin/env python3
# Writes the <pose>_score.csv of --report-format csv, in the current directory.
import sys
from pathlib import Path
pose = Path(sys.argv[sys.argv.index("-l") + 1])
assert Path(sys.argv[sys.argv.index("-p") + 1]).name.endswith("_prep.pdb"), "the receptor must be the charged PDB"
Path(pose.name.split(".")[0] + "_score.csv").write_text("ligand,DockTScore\\n" + pose.name + ",-6.5\\n")
"""


def executable(path: Path, text: str) -> Path:
    path.write_text(text)
    path.chmod(0o755)
    return path


def test_qualification_set_rescores_and_reports_the_scoring_power(tmp_path, capsys):
    """--docktscore-path: the score of the best-scored pose lands in the summary, next to the measured affinity of
    the complex, and scoring_power.txt is written beside success_rate.txt."""
    root = make_set(tmp_path)
    # A receptor the charged PDB can be written from, and the affinities a set keeps.
    (root / "1abc/receptors/1abc.in").write_text(
        "    1   34   -9.68400   13.42400   -1.74500      2     5     6     7     0     0 -0.853  N   GLY X   1    1\n"
    )
    (root / "index/affinity.csv").write_text("complex,log_affinity,delta_g_kcal_mol\n1abc,5.85,-7.979\n")
    grid_toml_generator.main([str(root)])
    args = [
        str(root),
        "-o",
        "runs",
        "--dockthor-path",
        str(executable(tmp_path / "dockthor", FAKE_DOCKTHOR)),
        "--docktscore-path",
        str(executable(tmp_path / "docktscore", FAKE_DOCKTSCORE)),
    ]
    capsys.readouterr()

    assert qualification_set.main(args) == 0
    row = read_summary(root / "runs" / "qualification_summary.csv")[0]
    assert row["complex"] == "1abc"
    assert [row[column_of(model)] for model in DEFAULT_MODELS] == ["-6.5"] * len(DEFAULT_MODELS)
    assert (row["log_affinity"], row["delta_g_kcal_mol"]) == ("5.85", "-7.979")
    # The pose scored is the first of the docked file, and its files stay with the complex.
    assert "best" in (root / "runs/1abc/1abc_top1.mol2").read_text()
    assert "second" not in (root / "runs/1abc/1abc_top1.mol2").read_text()
    assert (root / "runs/1abc/1abc_top1_general_dock_score.csv").is_file()
    assert (root / "runs/1abc/1abc_scores.csv").is_file()  # the one row a parallel job hands back
    # The charged PDB was written from the .in, into the set.
    assert (root / "1abc/receptors/prep/1abc_prep.pdb").is_file()
    power = (root / "runs" / "scoring_power.txt").read_text()
    assert "Complexes with results: 1" in power
    assert "Pearson r      : N/A" in power  # one complex: nothing to correlate, and it says so
    for model in DEFAULT_MODELS:  # one block per model, and a row per model in the comparison table
        assert f"## {model} ({column_of(model)})" in power
    scoring_csv = read_summary(root / "runs" / "scoring_power.csv")
    assert [row["model"] for row in scoring_csv] == DEFAULT_MODELS


def test_qualification_set_without_rescoring_writes_no_scoring_power(tmp_path, capsys):
    """Without --docktscore-path nothing is rescored, so there is no scoring power to report -- but the summary
    still carries the measured affinity."""
    root = make_set(tmp_path)
    (root / "index/affinity.csv").write_text("complex,log_affinity,delta_g_kcal_mol\n1abc,5.85,-7.979\n")
    grid_toml_generator.main([str(root)])
    capsys.readouterr()

    args = [str(root), "-o", "runs", "--dockthor-path", str(executable(tmp_path / "dockthor", FAKE_DOCKTHOR))]
    assert qualification_set.main(args) == 0
    row = read_summary(root / "runs" / "qualification_summary.csv")[0]
    assert not [column for column in row if column.startswith("docktscore")]
    assert row["delta_g_kcal_mol"] == "-7.979"
    assert not (root / "runs" / "scoring_power.txt").exists()
    with pytest.raises(SystemExit):  # --docktscore-args without the binary
        qualification_set.main([*args, "--docktscore-args=-v"])


def test_summarize_results(tmp_path, capsys, monkeypatch):
    reports = tmp_path / "reports"
    for pdb, top1 in (("1abc", "1.2"), ("2xyz", "3.4")):
        (reports / pdb).mkdir(parents=True)
        (reports / pdb / f"{pdb}_ligand_docked_success_rate.csv").write_text(
            f"complex,execution_time(s),num_top,top1_E,top1_RMSD\n{pdb},1,1,-9,{top1}\n"
        )
    rates = [
        "        0.5     0.00",
        "          1     0.00",
        "        1.5    50.00",
        "          2    50.00",
    ]
    assert summarize_results.main([str(reports)]) == 0
    assert (reports / "success_rate.txt").read_text().splitlines()[-4:] == rates
    assert read_summary(reports / "qualification_summary.csv")[1]["top1_RMSD"] == "3.4"

    # --config-dir counts against its TOMLs: 3def has no report.
    config = tmp_path / "config"
    config.mkdir()
    for pdb in ("1abc", "2xyz", "3def"):
        (config / f"{pdb}.toml").write_text("")
    assert summarize_results.main([str(reports), "--config-dir", str(config)]) == 0
    assert "1 complex(es) without results: 3def" in capsys.readouterr().out

    # --from-summary: only success_rate.txt, in the current directory.
    monkeypatch.chdir(tmp_path)
    assert summarize_results.main(["--from-summary", str(reports / "qualification_summary.csv")]) == 0
    assert (tmp_path / "success_rate.txt").read_text().splitlines()[-4:] == rates
    with pytest.raises(SystemExit):
        summarize_results.main(["--from-summary", "x.csv", "--config-dir", str(config)])
