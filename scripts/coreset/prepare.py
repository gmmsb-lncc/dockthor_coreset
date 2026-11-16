"""Build a qualification set from the CASF-2013 core set, with the DockThor Suite facilities.

The CASF-2013 distribution holds, per complex `<pdb>`:

    <pdb>/<pdb>_protein.pdb        the receptor, hydrogens included
    <pdb>/<pdb>_ligand.mol2        the ligand in the deposited (crystal) pose
    <pdb>/<pdb>_ligand.sdf         the same pose as SDF
    <pdb>/<pdb>_ligand_opt.mol2    the energy-minimized pose -- never used here (README 1.2)

and this module turns it into the layout of the set (README 1.4), one directory per complex: `<pdb>/receptors/<pdb>.in`
with PDBThorBox, `<pdb>/ligands/<pdb>_ligand.top` and `<pdb>/cofactors/` with MMFFLigand.

Two quirks of the CASF protein files are handled before PDBThorBox sees them, and both make it fail otherwise:

  * **NUL bytes in the chain column.** Several files (`1hfs`, `2jdu`, `2jdy`, ...) carry a 0x00 where the
    chain id belongs, which PDBThorBox reads as an extra, nameless chain.
  * **Metal ions written as `ATOM` records** whose residue name is the element (`ATOM ... ZN   ZN   1`).
    They have no residue template, so the chain cannot be built. They are the cofactors of the set, and
    `split_protein` moves them out of the receptor and into one file each.

Crystal waters are dropped, as they are in the committed set.
"""

import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from coreset import affinity
from coreset.dataset import COFACTORS_DIR, ligand_path, receptor_path

#: Residue names of the monatomic ions carried as cofactors. A record counts as one only when its atom name
#: is the residue name too, so the alpha carbon CA of an amino acid is never taken for a calcium.
ION_RESIDUES = frozenset("ZN CA MG NA MN HG CD NI FE CU CO K LI SR BA PB PT AU AG CS RB CL BR IOD F".split())
#: Dropped from the receptor, as in the committed set.
WATER_RESIDUES = frozenset(("HOH", "WAT", "DOD"))

STANDARD_RESIDUES = frozenset("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())


def is_ion(record: str) -> bool:
    """A PDB record holding one monatomic ion: residue name in ION_RESIDUES and the atom name the same."""
    residue = record[17:20].strip().upper()
    return residue in ION_RESIDUES and residue not in STANDARD_RESIDUES and record[12:16].strip().upper() == residue


@dataclass
class Split:
    """What `split_protein` found in one CASF protein file."""

    protein: list[str] = field(default_factory=list)
    ions: list[str] = field(default_factory=list)
    waters: int = 0
    nul_bytes: int = 0


def split_protein(pdb_file: Path) -> Split:
    """Separate a CASF protein file into the receptor records, the ion records and the dropped waters."""
    text = pdb_file.read_text(errors="replace")
    split = Split(nul_bytes=text.count("\0"))
    for line in text.replace("\0", " ").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            if line.startswith(("TER", "END")):
                continue
            split.protein.append(line)
            continue
        residue = line[17:20].strip().upper()
        if residue in WATER_RESIDUES:
            split.waters += 1
        elif is_ion(line):
            split.ions.append(line)
        else:
            split.protein.append(line)
    return split


@dataclass
class Tools:
    """The facility binaries used to prepare a complex; None for a facility that is skipped."""

    mmffligand: Path | None
    pdbthorbox: Path | None

    @staticmethod
    def find(bin_dir: Path, mmffligand: bool = True, pdbthorbox: bool = True) -> "Tools":
        """The binaries of `bin_dir`, only those asked for. Raise FileNotFoundError when one of them is missing."""
        tools = Tools(
            mmffligand=bin_dir / "mmffligand" if mmffligand else None,
            pdbthorbox=bin_dir / "pdbthorbox" if pdbthorbox else None,
        )
        missing = [str(p) for p in (tools.mmffligand, tools.pdbthorbox) if p is not None and not p.is_file()]
        if missing:
            raise FileNotFoundError("not a DockThor Suite bin directory, missing " + ", ".join(missing))
        return tools


class StepFailed(RuntimeError):
    """A facility exited non-zero; the message names the log file that holds its output."""


def run(command: list[str], log: Path, cwd: Path | None = None) -> None:
    """Run one facility, saving its console output to `log`. Raise StepFailed on a non-zero exit."""
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w") as stream:
        stream.write("$ " + " ".join(command) + "\n\n")
        stream.flush()
        status = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, cwd=cwd).returncode
    if status != 0:
        raise StepFailed(f"{Path(command[0]).name} exited {status}, see {log}")


@dataclass
class Prepared:
    """What one complex contributed to the set."""

    pdb: str
    ligand_atoms: int | None  # None: MMFFLigand skipped, the set's ligand and cofactors kept
    receptor_atoms: int | None  # None: PDBThorBox skipped, the set's receptor kept
    cofactors: int | None
    waters_dropped: int
    carried_cofactors: int = 0
    generated_files: int = 0


def count_top_atoms(top: Path) -> int:
    for line in top.read_text().splitlines():
        if line.startswith("$TOTAL_ATOMS"):
            return int(line.split()[1])
    raise ValueError(f"{top}: no $TOTAL_ATOMS")


def install(source: Path, destination: Path) -> None:
    """Copy one prepared file into the set, creating its folder."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def prepare_ligand(tools: Tools, casf: Path, pdb: str, out: Path, logs: Path) -> Path:
    """MMFFLigand on the CASF ligand, in the deposited pose. The `<pdb>_ligand.top` it wrote."""
    source = casf / pdb / f"{pdb}_ligand.mol2"
    if not source.is_file():
        raise StepFailed(f"missing {source}")
    # MMFFLigand writes beside its input, so the input is copied into the complex's output folder first.
    shutil.copy2(source, out / source.name)
    log = logs / f"{pdb}_mmffligand.log"
    run([str(tools.mmffligand), "-l", source.name], log, cwd=out)
    top = out / f"{pdb}_ligand.top"
    if not top.is_file():
        raise StepFailed(f"MMFFLigand wrote no {top.name}, see {log}")
    return top


def prepare_receptor(tools: Tools, split: Split, pdb: str, out: Path, logs: Path) -> Path:
    """PDBThorBox on the protein records of the split. The `<pdb>_protein.in` it wrote."""
    cleaned = out / f"{pdb}_protein.pdb"
    cleaned.write_text("\n".join(split.protein) + "\nEND\n")
    log = logs / f"{pdb}_pdbthorbox.log"
    # Everything PDBThorBox writes -- the topology, the charged PDB and MOL2, the chain files, the residue
    # and summary reports, resumo.out -- stays in the complex's output folder.
    run([str(tools.pdbthorbox), "-r", cleaned.name, "-o", "."], log, cwd=out)
    topology = out / f"{pdb}_protein.in"
    if not topology.is_file():
        raise StepFailed(f"PDBThorBox wrote no {topology.name}, see {log}")
    return topology


#: Replaces a coordinate of exactly 0 when MMFFLigand rejects an ion (see `nudge_zero_coordinates`).
NUDGE = 0.001


def nudge_zero_coordinates(record: str) -> str | None:
    """A copy of `record` with every coordinate of exactly 0 moved by NUDGE A, or None when it has none.

    MMFFLigand rejects a molecule for which OpenBabel's `OBMol::Has3D()` is false, and for a single atom that
    is true only when none of its three coordinates is exactly 0 -- so a lone ion on an axis or a plane of the
    crystal frame is refused as "not 3D structure". Three ions of the set are: the calcium of `1f8b`, `1f8c`
    and `1f8d`, at x = y = 0. Moving them by 0.001 A makes MMFFLigand accept them, and reproduces the
    coordinates of the committed cofactors, which carry the same workaround from their 2019 preparation.
    """
    values = [float(record[30:38]), float(record[38:46]), float(record[46:54])]
    if not any(value == 0.0 for value in values):
        return None
    moved = "".join(f"{(NUDGE if value == 0.0 else value):8.3f}" for value in values)
    return record[:30] + moved + record[54:]


def prepare_cofactors(tools: Tools, split: Split, pdb: str, out: Path, logs: Path) -> list[Path]:
    """One topology per ion of the split, in `<out>/cofactors/`. The `.top` files written, in order.

    They go in a subfolder of their own because MMFFLigand rewrites `mmffligand_report.csv` in its working
    directory on every call, and the ligand's report must survive.
    """
    if not split.ions:
        return []
    work = out / "cofactors"
    work.mkdir(parents=True, exist_ok=True)
    tops = []
    for index, record in enumerate(split.ions):
        name = f"{pdb}_cofactor-{index:02d}"
        pdb_file = work / f"{name}.pdb"
        top = work / f"{name}.top"
        log = logs / f"{name}.log"
        pdb_file.write_text(record.rstrip() + "\n")
        run([str(tools.mmffligand), "-l", pdb_file.name], log, cwd=work)
        if not top.is_file() and (nudged := nudge_zero_coordinates(record)) is not None:
            pdb_file.write_text(nudged.rstrip() + "\n")
            run([str(tools.mmffligand), "-l", pdb_file.name], log, cwd=work)
        if not top.is_file():
            raise StepFailed(f"MMFFLigand wrote no {top.name}, see {log}")
        tops.append(top)
    return tops


def carry_cofactors(source: Path, pdb: str, dataset: Path, start: int) -> int:
    """Copy the cofactors of `source/<pdb>/cofactors/` (another set) that the CASF protein file does not hold,
    numbered after the ones already written. Used for the ions of `1n1m`, `3ao4` and `3cj2` (README 1.6)."""
    folder = source / pdb / COFACTORS_DIR
    if not folder.is_dir():
        return 0
    tops = sorted(folder.glob("*.top"))
    for index, top in enumerate(tops, start=start):
        name = f"{pdb}_cofactor-{index:02d}"
        install(top, dataset / pdb / COFACTORS_DIR / f"{name}.top")
        companion = top.with_suffix(".pdb")
        if companion.is_file():
            install(companion, dataset / pdb / COFACTORS_DIR / f"{name}.pdb")
    return len(tops)


def prepare_complex(
    tools: Tools,
    casf: Path,
    pdb: str,
    output_dir: Path,
    dataset: Path,
    extra_cofactors: Path | None = None,
) -> Prepared:
    """Prepare one complex with the facilities of `tools` (a None one is skipped).

    `output_dir/<pdb>/` keeps **everything** the facilities wrote, one folder per complex, with the console
    output of every call in `output_dir/<pdb>/logs/`. `dataset` then receives only the files the set is made
    of (README 1.4), in `dataset/<pdb>/`: the ligand topology and the cofactor topologies (MMFFLigand) and the
    receptor topology under the name of the complex (PDBThorBox). What a skipped facility would have written is
    left in the set as it is. Raises StepFailed when a facility fails.
    """
    protein = casf / pdb / f"{pdb}_protein.pdb"
    if not protein.is_file():
        raise StepFailed(f"missing {protein}")
    split = split_protein(protein)

    out = output_dir / pdb
    logs = out / "logs"
    logs.mkdir(parents=True, exist_ok=True)

    top = None
    cofactor_tops: list[Path] = []
    if tools.mmffligand is not None:
        top = prepare_ligand(tools, casf, pdb, out, logs)
        cofactor_tops = prepare_cofactors(tools, split, pdb, out, logs)
    topology = prepare_receptor(tools, split, pdb, out, logs) if tools.pdbthorbox is not None else None

    if top is not None:
        install(top, dataset / ligand_path(pdb))
    if topology is not None:
        install(topology, dataset / receptor_path(pdb))
    for cofactor in cofactor_tops:
        install(cofactor, dataset / pdb / COFACTORS_DIR / cofactor.name)
        install(cofactor.with_suffix(".pdb"), dataset / pdb / COFACTORS_DIR / f"{cofactor.stem}.pdb")

    carried = 0
    if extra_cofactors is not None and tools.mmffligand is not None and not cofactor_tops:
        carried = carry_cofactors(extra_cofactors, pdb, dataset, start=len(cofactor_tops))

    return Prepared(
        pdb=pdb,
        ligand_atoms=count_top_atoms(top) if top is not None else None,
        receptor_atoms=(
            sum(1 for line in topology.read_text().splitlines() if line.split()) if topology is not None else None
        ),
        cofactors=len(cofactor_tops) if tools.mmffligand is not None else None,
        waters_dropped=split.waters,
        carried_cofactors=carried,
        generated_files=sum(1 for f in out.rglob("*") if f.is_file()),
    )


def write_index(casf: Path, dataset: Path) -> int:
    """Write `<dataset>/index/affinity.csv` from CASF's index files, and keep those files beside it.

    The index is the only part of the CASF distribution a set keeps: the affinity of each complex is measured
    data, not something the facilities can produce. CASF's two `.lst` are copied verbatim as the provenance of
    the numbers; the CSV converted from them is what every script reads (coreset/affinity.py).

    The number of complexes written, 0 when the CASF directory holds no index.
    """
    source = casf / "index"
    if not (source / affinity.CASF_DATA_FILE).is_file():
        return 0
    for name in (affinity.CASF_DATA_FILE, affinity.CASF_NAME_FILE):
        if (source / name).is_file():
            install(source / name, dataset / "index" / name)
    entries = affinity.convert_casf_index(source)
    affinity.write_csv(entries, dataset / "index" / affinity.CSV_FILE)
    return len(entries)


def discover_casf(casf: Path) -> list[str]:
    """Every complex of a CASF core-set directory: a folder holding both a protein and a ligand file."""
    return sorted(
        d.name
        for d in casf.iterdir()
        if d.is_dir() and (d / f"{d.name}_protein.pdb").is_file() and (d / f"{d.name}_ligand.mol2").is_file()
    )
