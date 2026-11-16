"""coreset.affinity: the CASF index converted to the `index/affinity.csv` a set keeps."""

import csv

import pytest

from coreset.affinity import (
    KCAL_PER_LOG_UNIT,
    NotAnIndex,
    columns,
    convert_casf_index,
    delta_g,
    parse_measurement,
    read_casf_names,
    write_csv,
)

DATA = """#===============================================================================
# The PDBbind core set v.2013, latest update: Nov 2013
# PDB code, resolution, release year, -logKd/Ki, original Kd/Ki, cluster ID
#===============================================================================
1ps3   1.80  2003   2.28  Ki=5.2mM     1
3d4z   1.39  2008   4.89  Ki=13uM      1
1bcu   2.00  1998   3.28  Kd=0.53mM   32
"""

NAMES = """#===============================================================================
# PDB code, release year, UniProt ID, protein name
#===============================================================================
1ps3  2003  Q24451   ALPHA-MANNOSIDASE II
1bcu  1998  P00734   THROMBIN
"""


def write_casf(root):
    """A CASF index directory holding both `.lst` files."""
    index = root / "index"
    index.mkdir(parents=True)
    (index / "2013_core_data.lst").write_text(DATA)
    (index / "2013_core_name.lst").write_text(NAMES)
    return index


def test_convert_reads_every_row_and_merges_the_protein_names(tmp_path):
    entries = {entry.complex: entry for entry in convert_casf_index(write_casf(tmp_path))}
    assert sorted(entries) == ["1bcu", "1ps3", "3d4z"]
    entry = entries["1bcu"]
    assert (entry.resolution_angstrom, entry.year, entry.log_affinity, entry.cluster) == (2.0, 1998, 3.28, 32)
    assert (entry.measurement, entry.value_molar, entry.original) == ("Kd", 0.00053, "Kd=0.53mM")
    assert entry.delta_g_kcal_mol == pytest.approx(-4.474, abs=1e-3)
    assert (entry.uniprot, entry.protein) == ("P00734", "THROMBIN")
    # A complex missing from the name file keeps its row, with no name.
    assert (entries["3d4z"].uniprot, entries["3d4z"].protein) == ("", "")


def test_delta_g_is_the_conversion_casf_uses():
    assert delta_g(3.28) == pytest.approx(-4.474, abs=1e-3)
    assert KCAL_PER_LOG_UNIT == -1.364


@pytest.mark.parametrize(
    ("text", "expected"),
    [("Ki=5.2mM", ("Ki", 0.0052)), ("Kd=0.53mM", ("Kd", 0.00053)), ("Ki=13uM", ("Ki", 1.3e-05))],
)
def test_parse_measurement(text, expected):
    assert parse_measurement(text) == expected


@pytest.mark.parametrize("text", ["5.2mM", "Ki=5.2", "Ki=5.2kM", "IC50=abcnM"])
def test_parse_measurement_rejects(text):
    with pytest.raises(NotAnIndex):
        parse_measurement(text)


def test_a_row_whose_constant_disagrees_with_its_log_is_rejected(tmp_path):
    index = write_casf(tmp_path)
    # 5.2 mM is -log10 2.28, not 8.28: a misread unit or column looks exactly like this.
    (index / "2013_core_data.lst").write_text(DATA.replace("2.28  Ki=5.2mM", "8.28  Ki=5.2mM"))
    with pytest.raises(NotAnIndex, match="8.28"):
        convert_casf_index(index)


def test_write_csv_one_row_per_complex(tmp_path):
    entries = convert_casf_index(write_casf(tmp_path))
    path = tmp_path / "set" / "index" / "affinity.csv"
    write_csv(entries, path)
    # Read back with the standard csv module, as a consumer of the set does.
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        assert reader.fieldnames == columns()
        rows = {row["complex"]: row for row in reader}
    assert sorted(rows) == sorted(entry.complex for entry in entries)
    assert rows["1bcu"]["protein"] == "THROMBIN"
    assert rows["3d4z"]["cluster"] == rows["1ps3"]["cluster"] == "1"  # one target, one cluster


def test_read_casf_names_skips_a_row_without_a_name(tmp_path):
    path = tmp_path / "n.lst"
    path.write_text(NAMES + "2xyz  2010  P12345\n")
    names = read_casf_names(path)
    assert "2xyz" not in names and names["1bcu"] == ("P00734", "THROMBIN")


@pytest.mark.parametrize(
    "text", ["# only a header\n", "1ps3 1.80 2003\n", "1ps3 x 2003 2.28 Ki=5mM 1\n", "1ps3 1.8 2003 2.28 5mM 1\n"]
)
def test_not_an_index(tmp_path, text):
    index = tmp_path / "index"
    index.mkdir()
    (index / "2013_core_data.lst").write_text(text)
    with pytest.raises(NotAnIndex):
        convert_casf_index(index)
