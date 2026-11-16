"""coreset.prepare: splitting a CASF protein file and the ion workaround."""

import pytest

from coreset.prepare import carry_cofactors, is_ion, nudge_zero_coordinates, split_protein

# A CASF-style protein file: a header, two residues, an ion as an ATOM record whose residue name is the
# element, an ion as a HETATM, a NUL byte in the chain column of one record, a water, and an alpha carbon
# that must not be taken for a calcium.
PROTEIN = (
    "HEADER    1ABC_PROTEIN\n"
    "ATOM      1  N   GLY A   1      11.000  12.000  13.000  1.00  0.00           N\n"
    "ATOM      2  CA  GLY A   1      11.500  12.500  13.500  1.00  0.00           C\n"
    "ATOM      3 ZN    ZN\x00   1      -4.317  -1.328  -4.742  1.00  0.00\n"
    "HETATM    4 CA    CA     2       3.100  11.613 -13.450  1.00  0.00\n"
    "HETATM    5  O   HOH B   3       1.000   2.000   3.000  1.00  0.00           O\n"
    "TER\n"
    "END\n"
)


def test_split_moves_the_ions_out_and_drops_the_water(tmp_path):
    path = tmp_path / "1abc_protein.pdb"
    path.write_text(PROTEIN)
    split = split_protein(path)
    assert split.nul_bytes == 1
    assert split.waters == 1
    assert [record[12:16].strip() for record in split.ions] == ["ZN", "CA"]
    # The protein keeps the header and the two amino-acid records, the alpha carbon among them.
    assert [line[:6] for line in split.protein] == ["HEADER", "ATOM  ", "ATOM  "]
    assert any(line[12:16].strip() == "CA" and line[17:20].strip() == "GLY" for line in split.protein)
    # The NUL byte is gone from the record that carried it.
    assert all("\0" not in record for record in split.ions)


@pytest.mark.parametrize(
    ("record", "expected"),
    [
        ("ATOM      3 ZN    ZN     1      -4.317  -1.328  -4.742", True),
        ("HETATM    4 CA    CA     2       3.100  11.613 -13.450", True),
        # An alpha carbon of an amino acid, and a delta carbon: the residue name decides.
        ("ATOM      2  CA  GLY A   1      11.500  12.500  13.500", False),
        ("ATOM      7  CD  LYS A   4      11.500  12.500  13.500", False),
        # A ligand atom named after its residue but not an ion.
        ("HETATM    8  C1  LIG A   5      11.500  12.500  13.500", False),
    ],
)
def test_is_ion(record, expected):
    assert is_ion(record) is expected


def test_nudge_only_touches_zero_coordinates():
    record = "HETATM23870 CA    CA     2       0.000   0.000  71.148  1.00  0.00"
    nudged = nudge_zero_coordinates(record)
    assert nudged is not None
    # The committed cofactor of 1f8b carries exactly these coordinates.
    assert nudged[30:54] == "   0.001   0.001  71.148"
    assert nudged[:30] == record[:30] and nudged[54:] == record[54:]


def test_nudge_returns_none_without_a_zero():
    assert nudge_zero_coordinates("HETATM 2409 ZN    ZN     1      -4.317  -1.328  -4.742  1.00  0.00") is None


def test_carry_cofactors_from_another_set(tmp_path):
    # The Hg of 1n1m, kept by an earlier set as 1n1m/cofactors/<name>.{top,pdb}, lands in the new set's
    # 1n1m/cofactors/ with the set's numbering, after the ones already written.
    old = tmp_path / "old" / "1n1m" / "cofactors"
    old.mkdir(parents=True)
    (old / "hg.top").write_text("$TOTAL_ATOMS 1\n")
    (old / "hg.pdb").write_text("HETATM\n")
    dataset = tmp_path / "new"
    assert carry_cofactors(tmp_path / "old", "1n1m", dataset, start=2) == 1
    assert sorted(p.name for p in (dataset / "1n1m" / "cofactors").iterdir()) == [
        "1n1m_cofactor-02.pdb",
        "1n1m_cofactor-02.top",
    ]
    assert carry_cofactors(tmp_path / "old", "3ao4", dataset, start=0) == 0
