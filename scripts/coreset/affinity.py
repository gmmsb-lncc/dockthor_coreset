"""The binding affinities of a qualification set: `index/affinity.csv`.

The affinities are the one thing the preparation cannot regenerate — the structures come out of the facilities,
but a measured affinity is data — so a set keeps them. CASF distributes them as two fixed-column text files with
`#` headers (`2013_core_data.lst`, `2013_core_name.lst`), which is awkward to read; a set stores them as one CSV
instead, written by `convert_casf_index` and read with the standard `csv` module:

    complex,resolution_angstrom,year,log_affinity,delta_g_kcal_mol,measurement,value_molar,original,cluster,uniprot,protein
    1bcu,2.0,1998,3.28,-4.474,Kd,0.00053,Kd=0.53mM,32,P00734,THROMBIN

| Column | Meaning |
| ------ | ------- |
| `complex` | PDB code, the name the set uses everywhere |
| `resolution_angstrom`, `year` | Of the crystal structure |
| `log_affinity` | CASF's `-logKd/Ki`, the value the scoring power correlates against |
| `delta_g_kcal_mol` | `KCAL_PER_LOG_UNIT * log_affinity`, to compare with a computed energy |
| `measurement` | `Ki` or `Kd` |
| `value_molar` | The published constant in mol/L (`-log10` of it is `log_affinity`, checked on conversion) |
| `original` | The published constant as CASF prints it, e.g. `Kd=0.53mM` |
| `cluster` | 1..65; the three complexes of one target share it -- the grouping of the ranking power |
| `uniprot`, `protein` | From CASF's name file; empty when it has no row for the complex |

Only the conversion lives here (`scripts/qualify_coreset.py` runs it, through `coreset.prepare.write_index`); a
reader of the CSV, for the correlations with the scores, needs nothing more than `csv.DictReader`.
"""

import csv
import math
import re
from dataclasses import asdict, dataclass, fields
from pathlib import Path

from coreset.dataset import INDEX_DIR

#: ΔG = -1.364 × pKd at 25 °C, in kcal/mol: the conversion CASF-2013 uses to compare with computed energies.
KCAL_PER_LOG_UNIT = -1.364

#: The file a set stores, and the CASF files it is converted from.
CSV_FILE = "affinity.csv"
CASF_DATA_FILE = "2013_core_data.lst"
CASF_NAME_FILE = "2013_core_name.lst"

CASF_DATA_FIELDS = 6  # code, resolution, year, -logKd/Ki, Kd/Ki, cluster
CASF_NAME_FIELDS = 4  # code, year, UniProt id, protein name (the rest of the line)

#: mol/L per unit of a published constant.
MOLAR_UNITS = {"M": 1.0, "mM": 1e-3, "uM": 1e-6, "nM": 1e-9, "pM": 1e-12, "fM": 1e-15}
MEASUREMENT = re.compile(r"(?P<kind>Ki|Kd|IC50)=(?P<value>[0-9.]+)(?P<unit>[a-zA-Z]+)")
#: How far -log10(value_molar) may sit from CASF's own -logKd/Ki before the row is rejected. The 195 rows of
#: CASF-2013 all agree to better than 0.05.
LOG_TOLERANCE = 0.05


class NotAnIndex(ValueError):
    """A file is not the affinity index it was read as."""


@dataclass(frozen=True)
class Entry:
    """The affinity of one complex, one row of `index/affinity.csv`."""

    complex: str
    resolution_angstrom: float
    year: int
    log_affinity: float
    delta_g_kcal_mol: float
    measurement: str
    value_molar: float
    original: str
    cluster: int
    uniprot: str = ""
    protein: str = ""


def read_affinity_csv(path: Path) -> dict[str, dict[str, str]]:
    """{complex: {column: value}} of an index/affinity.csv, the values as the text of the file; empty if absent."""
    if not path.is_file():
        return {}
    with open(path, newline="") as stream:
        return {row["complex"]: row for row in csv.DictReader(stream) if row.get("complex")}


def read_affinities(dataset: Path) -> dict[str, dict[str, str]]:
    """The affinities of a set, by complex: `<dataset>/index/affinity.csv` (empty when the set has no index)."""
    return read_affinity_csv(dataset / INDEX_DIR / CSV_FILE)


def columns() -> list[str]:
    return [field.name for field in fields(Entry)]


def delta_g(log_affinity: float) -> float:
    """The binding free energy in kcal/mol of a `-logKd/Ki`."""
    return KCAL_PER_LOG_UNIT * log_affinity


def parse_measurement(text: str) -> tuple[str, float]:
    """`'Kd=0.53mM'` -> `('Kd', 0.00053)`: the kind and the constant in mol/L."""
    match = MEASUREMENT.fullmatch(text)
    if match is None:
        raise NotAnIndex(f"'{text}' is not a published constant (Ki=, Kd= or IC50= with a unit)")
    unit = match.group("unit")
    if unit not in MOLAR_UNITS:
        raise NotAnIndex(f"'{text}': unknown unit '{unit}'")
    molar = float(match.group("value")) * MOLAR_UNITS[unit]
    # Scaling by the unit leaves binary noise (5.2 * 1e-3 is 0.005200000000000001); the published constants have
    # three significant figures at most, so six keep the value and write it as 0.0052 and 1.3e-05.
    return match.group("kind"), float(f"{molar:.6g}")


def read_casf_names(path: Path) -> dict[str, tuple[str, str]]:
    """The UniProt id and protein name of a CASF `*_core_name.lst`, by PDB code."""
    names = {}
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.split():
            continue
        row = line.split(maxsplit=CASF_NAME_FIELDS - 1)
        if len(row) < CASF_NAME_FIELDS:
            continue  # a row without a protein name
        names[row[0]] = (row[2], row[3].strip())
    return names


def convert_casf_index(index_dir: Path) -> list[Entry]:
    """Read CASF's index files from `index_dir` and return the entries, in file order.

    Every row is checked: `-log10(value_molar)` must agree with CASF's own `-logKd/Ki` to LOG_TOLERANCE, which
    catches a misread unit or column.
    """
    data = index_dir / CASF_DATA_FILE
    if not data.is_file():
        raise NotAnIndex(f"missing {data}")
    names = read_casf_names(index_dir / CASF_NAME_FILE) if (index_dir / CASF_NAME_FILE).is_file() else {}

    entries = []
    for number, line in enumerate(data.read_text().splitlines(), start=1):
        if line.startswith("#") or not line.split():
            continue
        row = line.split()
        if len(row) < CASF_DATA_FIELDS:
            raise NotAnIndex(f"{data}:{number}: expected {CASF_DATA_FIELDS} fields, found {len(row)}")
        try:
            log_affinity = float(row[3])
            kind, molar = parse_measurement(row[4])
            if abs(-math.log10(molar) - log_affinity) > LOG_TOLERANCE:
                raise NotAnIndex(
                    f"{data}:{number}: {row[4]} is -log10 {-math.log10(molar):.2f}, "
                    f"but the file says {log_affinity}"
                )
            uniprot, protein = names.get(row[0], ("", ""))
            entries.append(
                Entry(
                    complex=row[0],
                    resolution_angstrom=float(row[1]),
                    year=int(row[2]),
                    log_affinity=log_affinity,
                    delta_g_kcal_mol=round(delta_g(log_affinity), 3),
                    measurement=kind,
                    value_molar=molar,
                    original=row[4],
                    cluster=int(row[5]),
                    uniprot=uniprot,
                    protein=protein,
                )
            )
        except ValueError as error:
            raise NotAnIndex(f"{data}:{number}: {error}") from None
    if not entries:
        raise NotAnIndex(f"{data}: no entry")
    return entries


def write_csv(entries: list[Entry], path: Path) -> None:
    """Write `index/affinity.csv`, one row per complex."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns())
        writer.writeheader()
        for entry in entries:
            writer.writerow(asdict(entry))
