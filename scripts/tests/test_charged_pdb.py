"""coreset.charged_pdb: the charged PDB DockTScore reads, written from a DockThor receptor topology (`.in`)."""

import pytest

from coreset.charged_pdb import NotATopology, read_atoms, record, write_charged_pdb

# Three atom lines of a receptor topology, as PDBThorBox writes them: index, MMFF type, x, y, z, six bonded
# neighbors, charge, name, residue, chain, residue number, flag. The last is residue 1000, where the fixed columns
# of the `.in` leave no space between the chain and the number.
TOPOLOGY = (
    "    1   34   -9.68400   13.42400   -1.74500      2     5     6     7     0     0 -0.853  N   GLY X   1    1\n"
    "   49   28    1.00000    2.00000    3.00000     48     0     0     0     0     0  0.450  HH11 ARG X   5    1\n"
    " 9861   10   53.75800   42.51900  -18.53000   9859  9862  9863     0     0     0 -0.730  N   GLU A1000    1\n"
)

# The same three atoms as PDBThorBox's own `<receptor>_prep.pdb` has them.
CHARGED = (
    "ATOM      1  N   GLY X   1      -9.684  13.424  -1.745 -0.85 34.00           N\n"
    "ATOM     49 HH11 ARG X   5       1.000   2.000   3.000  0.45 28.00           H\n"
    "ATOM   9861  N   GLU A1000      53.758  42.519 -18.530 -0.73 10.00           N\n"
)


def test_record_is_the_pdb_of_pdbthorbox(tmp_path):
    """Every record is the one PDBThorBox writes: the fixed columns, 78 characters, charge in the occupancy column
    and MMFF type in the temperature factor."""
    topology = tmp_path / "3ehy.in"
    topology.write_text(TOPOLOGY)
    out = tmp_path / "prep" / "3ehy_prep.pdb"

    assert write_charged_pdb(topology, out) == 3
    assert out.read_text() == CHARGED
    assert {len(line) for line in out.read_text().splitlines()} == {78}


@pytest.fixture
def atoms(tmp_path):
    """The three atoms of TOPOLOGY, read back from a file as read_atoms takes a path."""
    topology = tmp_path / "x.in"
    topology.write_text(TOPOLOGY)
    return read_atoms(topology)


def test_columns_of_a_record(atoms):
    """The fields land in the columns the PDB format gives them (counted from 0 here, so 54:60 is columns 55-60)."""
    line = record(atoms[0])
    assert line[:6] == "ATOM  "
    assert line[6:11] == "    1"  # serial
    assert line[12:16] == " N  "  # atom name, left-justified in 14-16 when shorter than four characters
    assert line[17:20] == "GLY"  # residue name
    assert line[21] == "X"  # chain
    assert line[22:26] == "   1"  # residue number
    assert line[30:54] == "  -9.684  13.424  -1.745"
    assert line[54:60] == " -0.85"  # occupancy = partial charge
    assert line[60:66] == " 34.00"  # temperature factor = MMFF94S atom type
    assert line[76:78] == " N"  # element


def test_four_character_atom_name_starts_one_column_earlier(atoms):
    """A name of four characters fills columns 13-16, as the PDB format has it."""
    assert record(atoms[1])[12:16] == "HH11"


def test_chain_and_residue_number_run_together_past_residue_999(atoms):
    """The `.in` is written in fixed columns, so residue 1000 reads as `A1000`: one field, split back in two."""
    assert (atoms[2].chain, atoms[2].residue_number) == ("A", "1000")
    assert record(atoms[2])[21:26] == "A1000"


@pytest.mark.parametrize(
    "line",
    [
        "    1   34   -9.68400   13.42400\n",  # too few fields
        "    1   34   x   13.42400   -1.74500      2     5     6     7     0     0 -0.853  N   GLY X   1    1\n",
        "",  # no atom line at all
    ],
)
def test_a_file_that_is_not_a_topology_is_refused(tmp_path, line):
    topology = tmp_path / "bad.in"
    topology.write_text(line)
    with pytest.raises(NotATopology):
        write_charged_pdb(topology, tmp_path / "out.pdb")
