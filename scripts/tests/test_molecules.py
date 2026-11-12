"""coreset.molecules: the ligand topologies DockThor reads, and the other files of a ligand folder."""

import gzip

import pytest

from coreset.molecules import NotALigand, find_ligands, read_coordinates

TOP = """$TOTAL_ATOMS  3
   N     1    34     -0.85300     15.08800     10.79800     23.54700     7
   C     2     1      0.39700     15.01000      9.98700     24.79200     6
   H     3     5      0.00000     14.00000      9.00000     25.00000     1
$CON
1 2
"""

CIF = """data_LIG
#
_software.name                           mmffligand
_dts_block.command                       'mmffligand -l lig.mol2 --molecule-format cif'
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_ordinal
_chem_comp_atom.model_Cartn_x
_chem_comp_atom.model_Cartn_y
_chem_comp_atom.model_Cartn_z
LIG C1 C C 1 8.982 23.182 49.516
LIG "N 2" N N 2 9.67 23.649 48.475
LIG H3 H H 3 10.665 22.851 47.794
#
loop_
_mmffligand_atom.atom_id
_mmffligand_atom.mmff_type
_mmffligand_atom.partial_charge
C1 1 0.1
'N 2' 34 -0.2
H3 5 0.1
#
_mmffligand_ligand.note
;
a text field
over two lines
;
"""


def write(path, text):
    path.write_text(text)
    return path


def test_top(tmp_path):
    atoms = read_coordinates(write(tmp_path / "lig.top", TOP))
    assert atoms == [(15.088, 10.798, 23.547), (15.01, 9.987, 24.792), (14.0, 9.0, 25.0)]


def test_cif_and_gzipped_cif(tmp_path):
    expected = [(8.982, 23.182, 49.516), (9.67, 23.649, 48.475), (10.665, 22.851, 47.794)]
    assert read_coordinates(write(tmp_path / "lig_mmff.cif", CIF)) == expected
    with gzip.open(tmp_path / "lig_mmff.mmcif.gz", "wt") as stream:
        stream.write(CIF)
    assert read_coordinates(tmp_path / "lig_mmff.mmcif.gz") == expected


@pytest.mark.parametrize(
    ("name", "text"),
    [
        ("entry.cif", CIF.replace("_mmffligand_atom", "_other_atom")),  # a wwPDB entry: no MMFF94S types
        ("short.top", TOP.replace("$TOTAL_ATOMS  3", "$TOTAL_ATOMS  5")),
        ("bad.top", "not a topology\n"),
        ("empty.top", ""),
        ("lig.mol2", "@<TRIPOS>MOLECULE\n"),
    ],
)
def test_not_a_ligand(tmp_path, name, text):
    with pytest.raises(NotALigand):
        read_coordinates(write(tmp_path / name, text))


def test_find_ligands(tmp_path):
    write(tmp_path / "a_ligand.top", TOP)
    write(tmp_path / "a_ligand_mmff.cif", CIF)  # a second format of a: the .top wins, the .cif is not reported
    write(tmp_path / "a_ligand_mmff.mol2", "@<TRIPOS>MOLECULE\n")  # not read by DockThor
    write(tmp_path / "b_mmff.cif", CIF)
    write(tmp_path / "c.TOP", "garbage\n")  # a broken .top: the .cif of c is used instead
    write(tmp_path / "c_mmff.cif", CIF)
    write(tmp_path / "new_a_ligand.mol2", "x")  # deprecated copies: ignored without a message
    write(tmp_path / "a_ligand_new.mol2", "x")
    write(tmp_path / "a_ligand.top_diedral.inf", "x")
    (tmp_path / "subdir.top").mkdir()

    ligands, skipped = find_ligands(tmp_path)
    assert [(ligand.name, ligand.path.name) for ligand in ligands] == [
        ("a_ligand", "a_ligand.top"),
        ("b", "b_mmff.cif"),
        ("c", "c_mmff.cif"),
    ]
    assert [path.name for path in skipped] == ["a_ligand.top_diedral.inf", "a_ligand_mmff.mol2", "c.TOP"]
