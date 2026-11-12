"""coreset.tomlfile: the TOML literals and the <pdb>.toml layout."""

import argparse
import tomllib

import pytest

from coreset.tomlfile import parse_set_option, toml_text, toml_value


def test_toml_values():
    assert toml_value(True) == "true"
    assert toml_value(30) == "30"
    assert toml_value(0.25) == "0.25"
    assert toml_value([19.274, 19.274]) == "[19.274,19.274]"
    assert toml_value('a "quoted" path') == '"a \\"quoted\\" path"'


def test_toml_text_layout_and_comment():
    text = toml_text("3ag9", {"ligand": "ligand/3ag9_ligand.top", "grid-max-points": 1520875}, {"grid-max-points": "x"})
    assert text == (
        '## DockThor qualification set: 3ag9\n\nligand="ligand/3ag9_ligand.top"\n# x\ngrid-max-points=1520875\n'
    )
    assert tomllib.loads(text) == {"ligand": "ligand/3ag9_ligand.top", "grid-max-points": 1520875}


def test_set_option():
    assert parse_set_option("ga-nrun=30") == ("ga-nrun", 30)
    assert parse_set_option("--seed=7") == ("seed", 7)
    assert parse_set_option('ga-evolution-type="GENERATIONAL_OCL"') == ("ga-evolution-type", "GENERATIONAL_OCL")
    assert parse_set_option("output-dir=runs/a") == ("output-dir", "runs/a")  # not a TOML literal: kept as text
    with pytest.raises(argparse.ArgumentTypeError):
        parse_set_option("no-value")
