"""coreset.scoring: DockTScore run on the best-scored pose of a docked complex."""

import os
import stat
from argparse import Namespace
from pathlib import Path

from coreset.dataset import prep_path
from coreset.scoring import (
    DEFAULT_MODELS,
    SCORE_COLUMN,
    SCORES_CSV,
    cofactor_files,
    column_of,
    read_score,
    receptor_pdb,
    score_command,
    score_complex,
    split_best_pose,
)

TOPOLOGY = (
    "    1   34   -9.68400   13.42400   -1.74500      2     5     6     7     0     0 -0.853  N   GLY X   1    1\n"
)

POSES = """@<TRIPOS>MOLECULE
pose_one
@<TRIPOS>ATOM
      1 C      0.0000  0.0000  0.0000 C.3
@<TRIPOS>MOLECULE
pose_two
@<TRIPOS>ATOM
      1 C      1.0000  1.0000  1.0000 C.3
"""


def make_set(root: Path, pdb: str = "3ehy", cofactors: int = 2) -> Path:
    """A one-complex set: the receptor topology, the cofactor records, and no charged PDB yet."""
    (root / pdb / "receptors").mkdir(parents=True)
    (root / pdb / "receptors" / f"{pdb}.in").write_text(TOPOLOGY)
    (root / pdb / "cofactors").mkdir()
    for index in range(cofactors):
        (root / pdb / "cofactors" / f"{pdb}_cofactor-{index:02d}.pdb").write_text("ATOM\n")
        (root / pdb / "cofactors" / f"{pdb}_cofactor-{index:02d}.top").write_text("topology\n")
    return root


def fake_docktscore(path: Path, score: str = "-6.946457969891734", returncode: int = 0) -> Path:
    """A stand-in for the binary: it writes the <pose>_score.csv DockTScore would write, in the current directory."""
    path.write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "from pathlib import Path\n"
        "pose = Path(sys.argv[sys.argv.index('-l') + 1])\n"
        "name = pose.name.split('.')[0] + '_score.csv'\n"
        f"Path(name).write_text('ligand,{SCORE_COLUMN}\\n' + pose.name + ',{score}\\n')\n"
        f"sys.exit({returncode})\n"
    )
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return path


def test_split_best_pose_takes_the_first_of_the_poses(tmp_path):
    """DockThor writes the leaders best first, and DockTScore scores one pose at a time."""
    docked = tmp_path / "3ehy_ligand_docked.mol2"
    docked.write_text(POSES)
    pose = tmp_path / "3ehy_top1.mol2"

    assert split_best_pose(docked, pose)
    text = pose.read_text()
    assert text.count("@<TRIPOS>MOLECULE") == 1
    assert "pose_one" in text
    assert "pose_two" not in text


def test_split_best_pose_of_a_file_without_a_pose(tmp_path):
    empty = tmp_path / "empty.mol2"
    empty.write_text("")
    assert not split_best_pose(empty, tmp_path / "out.mol2")


def test_a_receptor_the_charged_pdb_cannot_be_written_from_is_a_warning(tmp_path, capsys):
    """A malformed .in must not stop the set: the complex is reported and left without a score."""
    dataset = make_set(tmp_path / "set")
    (dataset / "3ehy/receptors/3ehy.in").write_text("not a topology\n")
    out_dir = tmp_path / "out" / "3ehy"
    out_dir.mkdir(parents=True)
    (out_dir / "3ehy_ligand_docked.mol2").write_text(POSES)
    args = Namespace(
        dataset=dataset,
        docktscore_path=fake_docktscore(tmp_path / "docktscore.py"),
        docktscore_models=["general_dock"],
        docktscore_args=[],
    )

    assert score_complex(args, "3ehy", out_dir) == {}
    assert "no receptor to score against" in capsys.readouterr().err


def test_the_charged_pdb_is_written_from_the_in_when_the_set_has_none(tmp_path):
    """A set whose receptors/prep/ is missing still scores: the charged PDB is derived from the committed .in."""
    dataset = make_set(tmp_path / "set")
    charged = receptor_pdb(dataset, "3ehy")
    assert charged == dataset / prep_path("3ehy")
    assert charged.read_text().startswith("ATOM      1  N   GLY X   1")
    assert receptor_pdb(dataset, "3ehy") == charged  # kept, not written again


def test_only_the_pdb_cofactors_are_passed(tmp_path):
    """DockTScore reads a cofactor as PDB or MOL2, never as a `.top` topology."""
    dataset = make_set(tmp_path / "set")
    names = [path.name for path in cofactor_files(dataset, "3ehy")]
    assert names == ["3ehy_cofactor-00.pdb", "3ehy_cofactor-01.pdb"]


def test_column_of_a_model():
    assert column_of("general_dock") == "docktscore_general_dock"


def test_score_command(tmp_path):
    command = score_command(
        Path("/bin/docktscore"), Path("r.pdb"), Path("p.mol2"), [Path("c1.pdb"), Path("c2.pdb")], "ppi", ["-v"]
    )
    assert command[command.index("--model") + 1] == "ppi"
    assert command[:3] == ["/bin/docktscore", "-p", "r.pdb"]
    assert command[3:5] == ["-l", "p.mol2"]
    assert "--report-format" in command and "csv" in command
    assert command[command.index("--cofactor") + 1 : command.index("--cofactor") + 3] == ["c1.pdb", "c2.pdb"]
    assert command[-1] == "-v"
    assert "--cofactor" not in score_command(Path("d"), Path("r"), Path("p"), [], "ppi", [])


def test_read_score(tmp_path):
    path = tmp_path / "x_score.csv"
    path.write_text(f"ligand,{SCORE_COLUMN}\nx.mol2,-6.5\n")
    assert read_score(path) == "-6.5"  # the text of the file, so the full precision is kept
    assert read_score(tmp_path / "absent.csv") == ""


def test_score_complex_scores_the_best_pose_with_every_model(tmp_path, capsys):
    dataset = make_set(tmp_path / "set")
    out_dir = tmp_path / "out" / "3ehy"
    out_dir.mkdir(parents=True)
    (out_dir / "3ehy_ligand_docked.mol2").write_text(POSES)
    args = Namespace(
        dataset=dataset,
        docktscore_path=fake_docktscore(tmp_path / "docktscore.py"),
        docktscore_models=DEFAULT_MODELS,
        docktscore_args=[],
    )

    scores = score_complex(args, "3ehy", out_dir)
    assert scores == {column_of(model): "-6.946457969891734" for model in DEFAULT_MODELS}
    assert (out_dir / "3ehy_top1.mol2").is_file()
    for model in DEFAULT_MODELS:
        # Each model scores its own copy of the pose: DockTScore names its output after the ligand.
        assert (out_dir / f"3ehy_top1_{model}.mol2").is_file()
        assert (out_dir / f"3ehy_top1_{model}_score.csv").is_file()
    row = (out_dir / f"3ehy{SCORES_CSV}").read_text().splitlines()
    assert row[0] == "complex," + ",".join(column_of(model) for model in DEFAULT_MODELS)
    assert row[1].startswith("3ehy,-6.946457969891734")


def test_score_complex_warns_and_gives_no_score_when_it_cannot_run(tmp_path, capsys):
    dataset = make_set(tmp_path / "set")
    out_dir = tmp_path / "out" / "3ehy"
    out_dir.mkdir(parents=True)
    args = Namespace(
        dataset=dataset,
        docktscore_path=fake_docktscore(tmp_path / "docktscore.py", returncode=3),
        docktscore_models=["general_dock"],
        docktscore_args=[],
    )

    assert score_complex(args, "3ehy", out_dir) == {}  # no *_docked.mol2 at all
    assert "no *_docked.mol2" in capsys.readouterr().err

    (out_dir / "3ehy_ligand_docked.mol2").write_text(POSES)
    assert score_complex(args, "3ehy", out_dir) == {}  # the binary fails
    assert "docktscore --model general_dock failed (exit 3)" in capsys.readouterr().err
    assert not (out_dir / f"3ehy{SCORES_CSV}").exists()  # nothing to hand back


def test_the_fake_binary_is_executable(tmp_path):
    """The stand-in must really run, or the two tests above would pass for the wrong reason."""
    path = fake_docktscore(tmp_path / "docktscore.py")
    assert os.access(path, os.X_OK)
