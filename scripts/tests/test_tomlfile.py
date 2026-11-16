"""coreset.tomlfile: the TOML literals and the <pdb>.toml layout."""

import tomllib

from coreset.tomlfile import toml_text, toml_value


def test_toml_values():
    assert toml_value(True) == "true"
    assert toml_value(30) == "30"
    assert toml_value(0.25) == "0.25"
    assert toml_value([19.274, 19.274]) == "[19.274,19.274]"
    assert toml_value('a "quoted" path') == '"a \\"quoted\\" path"'


def test_toml_text_layout_and_comment():
    text = toml_text(
        "3ag9", {"ligand": "3ag9/ligands/3ag9_ligand.top", "grid-max-points": 1520875}, {"grid-max-points": "x"}
    )
    assert text == (
        '## DockThor qualification set: 3ag9\n\nligand="3ag9/ligands/3ag9_ligand.top"\n# x\ngrid-max-points=1520875\n'
    )
    assert tomllib.loads(text) == {"ligand": "3ag9/ligands/3ag9_ligand.top", "grid-max-points": 1520875}
