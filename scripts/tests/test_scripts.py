"""The two command lines end to end on a small set: grid_toml_generator.py, and qualification_set.py --generate and
--collect."""

import tomllib

import grid_toml_generator
import qualification_set
from test_molecules import CIF, TOP


def make_set(root):
    """Two complexes: 1abc with a .top and a cofactor, 2xyz with only MMFFLigand's .cif and no receptor."""
    for directory in ("ligand", "receptor", "cofactor/1abc", "config"):
        (root / directory).mkdir(parents=True)
    (root / "ligand/1abc_ligand.top").write_text(TOP)
    (root / "ligand/1abc_ligand_mmff.mol2").write_text("@<TRIPOS>MOLECULE\n")
    (root / "ligand/2xyz_ligand_mmff.cif").write_text(CIF)
    (root / "receptor/1abc.in").write_text("receptor\n")
    (root / "cofactor/1abc/1abc_cofactor-01.top").write_text(TOP)
    return root


def test_grid_toml_generator(tmp_path, capsys):
    root = make_set(tmp_path)
    assert grid_toml_generator.main([str(root), "--padding", "20"]) == 1  # 2xyz has no receptor
    out, err = capsys.readouterr()
    assert "Skipped unrecognized ligand/1abc_ligand_mmff.mol2" in err
    assert "2xyz: missing" in err
    assert not (root / "config/2xyz.toml").exists()

    text = (root / "config/1abc.toml").read_text()
    config = tomllib.loads(text)
    # extent 1.088 x 1.798 x 1.453 + 20: a cube of 21.798 Å, 88 points per axis
    assert config == {
        "output-dir": "vcoreset_run",
        "ligand": "ligand/1abc_ligand.top",
        "receptor": "receptor/1abc.in",
        "cofactors-path": "cofactor/1abc",
        "grid-size": [21.798, 21.798, 21.798],
        "grid-center": [14.544, 9.899, 24.2735],
    }  # no grid-rstep: --rstep not given, DockThor's default 0.25 counts the points
    assert text.startswith("## DockThor qualification set: 1abc\n\n")
    assert "681472" in out  # 88^3

    # An existing TOML is kept; --box-shape ligand and a large padding: over the limit, grid-max-points written.
    assert grid_toml_generator.main([str(root), "--systems", "1abc", "--padding", "40", "--box-shape", "ligand"]) == 0
    assert "exists (use --overwrite)" in capsys.readouterr().out
    assert (
        grid_toml_generator.main(
            [str(root), "--systems", "1abc", "--padding", "40", "--box-shape", "ligand", "--overwrite"]
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
        grid_toml_generator.main([str(root), "--systems", "1abc", "--padding", "20", "--rstep", "0.25", "--overwrite"])
        == 0
    )
    config = tomllib.loads((root / "config/1abc.toml").read_text())
    assert config["grid-rstep"] == 0.25
    assert "grid-max-points" not in config
    assert (
        grid_toml_generator.main([str(root), "--systems", "1abc", "--padding", "20", "--rstep", "0.2", "--overwrite"])
        == 0
    )
    config = tomllib.loads((root / "config/1abc.toml").read_text())
    assert config["grid-rstep"] == 0.2
    assert config["grid-max-points"] == 109**3  # 21.798 / 0.2 = 108.99 -> 109 points per axis


def test_qualification_set_generate_and_collect(tmp_path, capsys):
    root = make_set(tmp_path)
    (root / "config/1abc.grid_center.txt").write_text("1.5 2 3\n")
    (root / "config/1abc.grid_size.txt").write_text("19.274 19.274 19.274\n")
    (root / "config/1abc.rstep.txt").write_text("0.25\n")
    assert qualification_set.main([str(root), "--generate", "-o", "runs", "--set", "ga-nrun=30"]) == 0
    assert (root / "1abc.toml").read_text() == (
        "## DockThor qualification set: 1abc\n\n"
        'output-dir="runs"\nligand="ligand/1abc_ligand.top"\nreceptor="receptor/1abc.in"\n'
        'cofactors-path="cofactor/1abc"\ngrid-size=[19.274,19.274,19.274]\ngrid-center=[1.5,2.0,3.0]\n'
        "grid-rstep=0.25\nga-nrun=30\n"
    )

    reports = tmp_path / "reports"
    for pdb, top1 in (("1abc", "1.2"), ("2xyz", "3.4")):
        (reports / pdb).mkdir(parents=True)
        (reports / pdb / f"{pdb}_ligand_docked_success_rate.csv").write_text(
            f"complex,execution_time(s),num_top,top1_E,top1_RMSD\n{pdb},1,1,-9,{top1}\n"
        )
    assert qualification_set.main(["--collect", str(reports)]) == 0
    assert (reports / "success_rate.txt").read_text().splitlines()[-4:] == [
        "        0.5     0.00",
        "          1     0.00",
        "        1.5    50.00",
        "          2    50.00",
    ]
    assert qualification_set.read_summary(reports / "qualification_summary.csv")[1]["top1_RMSD"] == "3.4"
