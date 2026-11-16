"""The charged PDB DockTScore reads as its receptor, written from a DockThor receptor topology (`.in`).

DockTScore's `-p` takes a PDB whose **occupancy column is the partial charge and whose temperature factor is the
MMFF94S atom type** (`docktscore --help`); it refuses an `.in` topology. PDBThorBox writes that file beside the
topology (`<receptor>_prep.pdb`), but a set keeps only the `.in`, and the two must agree: DockThor docks into the
`.in`, so the receptor DockTScore scores has to be the same atoms, charges and types. This module therefore derives
the charged PDB **from the committed `.in` itself**, which holds everything a record needs:

    1504   32  -18.92100   13.93800  -12.22700   1502  0 0 0 0 0  -0.900  OXT TER X 158    1
    index  type  x           y          z        <- bonded ->     charge  name res chain number

becomes, in the fixed columns of a PDB ATOM record (78 characters, as PDBThorBox writes them):

    ATOM   1504  OXT TER X 158     -18.921  13.938 -12.227 -0.90 32.00           O

| Columns | Field | From the `.in` |
| ------- | ----- | -------------- |
| 7-11 | atom serial | the atom index |
| 13-16 | atom name | the name, in columns 13-16 when it has four characters, else left-justified in 14-16 |
| 18-20 | residue name | the residue name (PDBThorBox writes `TER` for a chain-terminating `OXT`) |
| 22 | chain | the chain id |
| 23-26 | residue number | the residue number |
| 31-54 | x, y, z | the coordinates, `%8.3f` |
| 55-60 | **occupancy = partial charge**, `%6.2f` | the charge |
| 61-66 | **temperature factor = MMFF94S atom type**, `%6.2f` | the type |
| 77-78 | element | the first letter of the atom name |

The two decimals of the charge are PDBThorBox's own rounding, kept so the file is byte for byte the one DockTScore
is validated against (`write_charged_pdb` reproduces every `_prep.pdb` of a PDBThorBox run; see
`scripts/tests/test_charged_pdb.py`).
"""

from dataclasses import dataclass
from pathlib import Path

#: Fields of an atom line of a `.in`: index, type, x, y, z, six bonded neighbors, charge, name, residue, chain,
#: residue number, and a last flag. The `.in` is written in fixed columns with one column for the chain and four for
#: the residue number, so from residue 1000 on the two run together (`A1000`) and the line has one field less; six
#: receptors of the CASF-2013 core set are numbered that far (README 4.5).
IN_FIELDS = 17
IN_FIELDS_MERGED = 16
NAME_COLUMNS = 4  # an atom name of four characters starts one column earlier, as the PDB format has it


class NotATopology(ValueError):
    """The file is not a DockThor receptor topology (`.in`)."""


@dataclass(frozen=True)
class Atom:
    """One atom of a receptor topology, the fields a charged PDB record needs."""

    index: int
    atom_type: int
    x: float
    y: float
    z: float
    charge: float
    name: str
    residue: str
    chain: str
    residue_number: str


def read_atoms(topology_in: Path) -> list[Atom]:
    """The atoms of a `.in`. Raise NotATopology when a line is not an atom line of one."""
    atoms = []
    for number, line in enumerate(topology_in.read_text().splitlines(), start=1):
        fields = line.split()
        if not fields:
            continue
        if len(fields) == IN_FIELDS:
            chain, residue_number = fields[14], fields[15]
        elif len(fields) == IN_FIELDS_MERGED and fields[14][1:].isdigit():
            chain, residue_number = fields[14][0], fields[14][1:]  # `A1000`: residue 1000 and beyond
        else:
            raise NotATopology(f"{topology_in}:{number}: {len(fields)} field(s), not the {IN_FIELDS} of an atom line")
        try:
            atoms.append(
                Atom(
                    index=int(fields[0]),
                    atom_type=int(fields[1]),
                    x=float(fields[2]),
                    y=float(fields[3]),
                    z=float(fields[4]),
                    charge=float(fields[11]),
                    name=fields[12],
                    residue=fields[13],
                    chain=chain,
                    residue_number=residue_number,
                )
            )
        except ValueError as error:
            raise NotATopology(f"{topology_in}:{number}: {error}") from None
    if not atoms:
        raise NotATopology(f"{topology_in}: no atom line")
    return atoms


def record(atom: Atom) -> str:
    """One ATOM record of the charged PDB, 78 characters."""
    name = f"{atom.name:<4}" if len(atom.name) >= NAME_COLUMNS else f" {atom.name:<3}"
    element = next((character for character in atom.name if character.isalpha()), "")
    return (
        f"ATOM  {atom.index:>5} {name} {atom.residue:>3} {atom.chain:1}{atom.residue_number:>4}    "
        f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}{atom.charge:6.2f}{float(atom.atom_type):6.2f}"
        f"{element:>12}"
    )


def write_charged_pdb(topology_in: Path, out: Path) -> int:
    """Write the charged PDB of a `.in` topology; the number of atoms written."""
    atoms = read_atoms(topology_in)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(record(atom) for atom in atoms) + "\n")
    return len(atoms)
