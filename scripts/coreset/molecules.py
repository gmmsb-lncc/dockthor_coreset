"""The ligand topologies DockThor reads, and their atom coordinates.

DockThor reads a ligand (dockthor/src/Topology_handler.cpp, read_ligand) as:

- a `.top` topology (MMFFLigand's <name>.top): `$TOTAL_ATOMS N`, then one line per atom,
  `element index mmff_type charge x y z atomic_number`;
- an mmCIF ligand topology (MMFFLigand's <name>_mmff.cif, also `.mmcif` and gzipped): a data block with
  `_mmffligand_atom` (the MMFF94S types and charges; DockThor refuses a block without it) and the coordinates in
  `_chem_comp_atom.model_Cartn_x/y/z`.

Any other file (a `.mol2`, a wwPDB entry, notes) is not a ligand topology. MMFFLigand's deprecated protonated copies
(`new_<name>`, `<name>_new.<ext>`) are not either, and are ignored without a message.
"""

import gzip
import re
from dataclasses import dataclass
from pathlib import Path

# The ligand formats of DockThor, in the order one is chosen when a ligand has several.
LIGAND_SUFFIXES = (".top", ".cif", ".mmcif", ".cif.gz", ".mmcif.gz")
MMFF_SUFFIX = "_mmff"  # MMFFLigand's <name>_mmff.cif (and .mol2)
CIF_COORDINATES = tuple(f"_chem_comp_atom.model_Cartn_{axis}" for axis in "xyz")
CIF_TYPES_CATEGORY = "_mmffligand_atom."


class NotALigand(ValueError):
    """The file is not a ligand topology DockThor reads."""


@dataclass(frozen=True)
class Ligand:
    name: str  # the file name without the format suffix and MMFFLigand's _mmff
    path: Path
    coordinates: list[tuple[float, float, float]]


def ligand_suffix(path: Path) -> str | None:
    """The DockThor ligand suffix of a file name (any case), longest first so .cif.gz is not taken as .gz."""
    name = path.name.lower()
    return next((s for s in sorted(LIGAND_SUFFIXES, key=len, reverse=True) if name.endswith(s)), None)


def ligand_name(path: Path) -> str:
    """lig.top -> lig, lig_mmff.cif.gz -> lig."""
    suffix = ligand_suffix(path)
    stem = path.name[: -len(suffix)] if suffix else path.stem
    return stem.removesuffix(MMFF_SUFFIX)


def is_deprecated_copy(path: Path) -> bool:
    """MMFFLigand's protonated copies of an input: new_<name> and <name>_new.<ext>."""
    return path.name.startswith("new_") or re.search(r"_new\.[^/]+$", path.name) is not None


def read_text(path: Path) -> str:
    if path.name.lower().endswith(".gz"):
        with gzip.open(path, "rt") as stream:
            return stream.read()
    return path.read_text()


def coordinate(text: str, path: Path, where: str) -> float:
    try:
        return float(text)
    except ValueError:
        raise NotALigand(f"{path}: {where}: '{text}' is not a coordinate") from None


def read_top_coordinates(path: Path) -> list[tuple[float, float, float]]:
    """The atoms of a .top: the `$TOTAL_ATOMS N` line, then N atom lines with x y z in columns 5-7."""
    lines = [line for line in read_text(path).splitlines() if line.strip()]
    header = lines[0].split() if lines else []
    if len(header) != 2 or header[0] != "$TOTAL_ATOMS" or not header[1].isdigit() or int(header[1]) == 0:
        raise NotALigand(f"{path}: no '$TOTAL_ATOMS N' first line")
    count = int(header[1])
    atoms = lines[1 : 1 + count]
    coordinates = []
    for number, line in enumerate(atoms, start=2):
        fields = line.split()
        if len(fields) < 7 or fields[0].startswith("$"):
            raise NotALigand(f"{path}: {count} atoms announced, atom line {number - 1} is '{line.strip()}'")
        coordinates.append(tuple(coordinate(v, path, f"atom {number - 1}") for v in fields[4:7]))
    if len(coordinates) != count:
        raise NotALigand(f"{path}: {count} atoms announced, {len(coordinates)} found")
    return coordinates


CIF_TOKEN = re.compile(r"""'.*?'(?=\s|$)|".*?"(?=\s|$)|\S+""")


def cif_tokens(text: str) -> list[str]:
    """The tokens of a CIF data block: quoted strings unquoted, `;` text fields joined, comments dropped."""
    tokens = []
    lines = iter(text.splitlines())
    for line in lines:
        if line.startswith(";"):
            field = [line[1:]]
            for more in lines:
                if more.startswith(";"):
                    break
                field.append(more)
            tokens.append("\n".join(field))
            continue
        for token in CIF_TOKEN.findall(line):
            if token.startswith("#"):
                break
            tokens.append(token[1:-1] if len(token) >= 2 and token[0] == token[-1] and token[0] in "'\"" else token)
    return tokens


def cif_items(text: str) -> dict[str, list[str]]:
    """{_category.item: values} of the first data block of an mmCIF file (one value for a single item)."""
    tokens = cif_tokens(text)
    items: dict[str, list[str]] = {}
    index = 0
    started = False
    while index < len(tokens):
        token = tokens[index]
        lower = token.lower()
        if lower.startswith("data_"):
            if started:
                break  # only the first block
            started = True
            index += 1
        elif lower == "loop_":
            index += 1
            tags = []
            while index < len(tokens) and tokens[index].startswith("_"):
                tags.append(tokens[index])
                index += 1
            values = []
            while (
                index < len(tokens)
                and not (tokens[index].startswith("_") or tokens[index].lower() in ("loop_", "stop_"))
                and not tokens[index].lower().startswith(("data_", "save_"))
            ):
                values.append(tokens[index])
                index += 1
            for column, tag in enumerate(tags):
                items[tag] = values[column :: len(tags)] if tags else []
        elif token.startswith("_") and index + 1 < len(tokens):
            items[token] = [tokens[index + 1]]
            index += 2
        else:
            index += 1
    return items


def read_cif_coordinates(path: Path) -> list[tuple[float, float, float]]:
    """The atoms of MMFFLigand's mmCIF ligand topology: _chem_comp_atom.model_Cartn_x/y/z, with _mmffligand_atom."""
    try:
        items = cif_items(read_text(path))
    except (OSError, UnicodeDecodeError, EOFError) as error:
        raise NotALigand(f"{path}: {error}") from None
    if not any(tag.startswith(CIF_TYPES_CATEGORY) for tag in items):
        raise NotALigand(f"{path}: no MMFF94S types ({CIF_TYPES_CATEGORY.rstrip('.')}), not a MMFFLigand topology")
    if any(tag not in items for tag in CIF_COORDINATES):
        raise NotALigand(f"{path}: no {CIF_COORDINATES[0].rsplit('_', 1)[0]}_x/y/z coordinates")
    columns = [items[tag] for tag in CIF_COORDINATES]
    if not columns[0] or len({len(column) for column in columns}) != 1:
        raise NotALigand(f"{path}: no atom coordinates")
    return [
        tuple(coordinate(v, path, f"atom {n}") for v in xyz)
        for n, xyz in enumerate(zip(*columns, strict=True), start=1)
    ]


def read_coordinates(path: Path) -> list[tuple[float, float, float]]:
    """The atom coordinates of a ligand topology; NotALigand when DockThor would not read it as one."""
    suffix = ligand_suffix(path)
    if suffix is None:
        raise NotALigand(f"{path}: not a DockThor ligand format ({', '.join(LIGAND_SUFFIXES)})")
    try:
        if suffix == ".top":
            return read_top_coordinates(path)
        return read_cif_coordinates(path)
    except (OSError, UnicodeDecodeError) as error:
        raise NotALigand(f"{path}: {error}") from None


def find_ligands(directory: Path) -> tuple[list[Ligand], list[Path]]:
    """The ligands of a directory (not its subdirectories), by name, and the files that are not ligand topologies.

    A ligand found in several formats is read from the first one of LIGAND_SUFFIXES that DockThor reads; the other
    copies are neither used nor reported. Deprecated copies (new_*, *_new.*) are ignored."""
    candidates: dict[str, list[Path]] = {}
    skipped = []
    for path in sorted(p for p in directory.iterdir() if p.is_file()):
        if is_deprecated_copy(path):
            continue
        if ligand_suffix(path) is None:
            skipped.append(path)
            continue
        candidates.setdefault(ligand_name(path), []).append(path)

    ligands = []
    for name, paths in sorted(candidates.items()):
        paths.sort(key=lambda p: LIGAND_SUFFIXES.index(ligand_suffix(p)))
        for path in paths:
            try:
                ligands.append(Ligand(name, path, read_coordinates(path)))
                break
            except NotALigand:
                skipped.append(path)
    return ligands, sorted(skipped)
