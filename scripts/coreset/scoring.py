"""Rescoring the docked poses with DockTScore, right after DockThor docked them or on a finished run's folder.

DockThor ranks the poses it writes with its own energy; DockTScore is the suite's scoring function, and its value is
what the scoring power of a set is measured with (`coreset/statistics.py`). For every complex that docked, this
module scores **the best-scored pose**, the first of `<ligand>_docked.mol2`, which is the pose CASF's scoring power
uses, once per DockTScore model:

    <output-dir>/<pdb>/<pdb>_top1.mol2               the pose alone, split out of the multi-pose file
    <output-dir>/<pdb>/<pdb>_top1_<model>.mol2       the pose as given to one model (DockTScore names its
    <output-dir>/<pdb>/<pdb>_top1_<model>.score      output after the ligand, so each model needs its own name)
    <output-dir>/<pdb>/<pdb>_top1_<model>_score.csv  the same terms at full precision (--report-format csv)
    <output-dir>/<pdb>/<pdb>_scores.csv              one row: the complex and its score per model

The last file is what a parallel job hands back (`scripts/ci/dock.sh`): it is small, holds no coordinates, and
`summarize_results.py` merges it into the summary like a success-rate report.

DockTScore needs its receptor as a charged PDB, which the set keeps in `<pdb>/receptors/prep/` and which is written
from the `.in` when it is missing (`coreset/charged_pdb.py`), and reads cofactors only as PDB or MOL2, which is what
`<pdb>/cofactors/*.pdb` holds.
"""

import csv
import shutil
import subprocess
import sys
from pathlib import Path

from coreset.charged_pdb import NotATopology, write_charged_pdb
from coreset.dataset import COFACTORS_DIR, prep_path, receptor_path

#: The line that starts a molecule in a MOL2 file: the poses of <ligand>_docked.mol2 are one block each, the first
#: of them the best-scored one.
MOL2_MOLECULE = "@<TRIPOS>MOLECULE"
#: The DockTScore models a qualification scores with (`docktscore --model`): the one fitted for docking poses, the
#: general set of 2764 complexes, and the two target families of the core set. Every one of them scores a pose in
#: about 0.1 s -- none needs an MSMS surface.
DEFAULT_MODELS = ["general_dock", "general2764", "protease", "ppi"]
#: The column of a <pose>_score.csv holding the score, in kcal/mol.
SCORE_COLUMN = "DockTScore"
#: Prefix of the summary column of a model: docktscore_general_dock, ... (`coreset/summary.py`).
COLUMN_PREFIX = "docktscore"
#: One row per complex, the scores of every model: what a parallel job hands back.
SCORES_CSV = "_scores.csv"


def column_of(model: str) -> str:
    """The summary column of a model: `general_dock` -> `docktscore_general_dock`."""
    return f"{COLUMN_PREFIX}_{model}"


def split_best_pose(docked_mol2: Path, out: Path) -> bool:
    """Write the first pose of a multi-pose MOL2 to `out`; False when the file holds no pose."""
    text = docked_mol2.read_text()
    blocks = [block for block in text.split(MOL2_MOLECULE) if block.strip()]
    if not blocks:
        return False
    out.write_text(MOL2_MOLECULE + blocks[0].rstrip("\n") + "\n")
    return True


def receptor_pdb(dataset: Path, pdb: str) -> Path:
    """The charged PDB of a complex, written from the committed `.in` when the set has none."""
    charged = dataset / prep_path(pdb)
    if not charged.is_file():
        write_charged_pdb(dataset / receptor_path(pdb), charged)
    return charged


def cofactor_files(dataset: Path, pdb: str) -> list[Path]:
    """The cofactor records of a complex, in name order (DockTScore reads a cofactor as PDB or MOL2, never `.top`)."""
    return sorted((dataset / pdb / COFACTORS_DIR).glob("*.pdb"))


def read_score(csv_path: Path) -> str:
    """The DockTScore value of a <pose>_score.csv, as the full-precision text of the file; empty if it has none."""
    if not csv_path.is_file():
        return ""
    with open(csv_path, newline="") as stream:
        rows = list(csv.DictReader(stream))
    return rows[0].get(SCORE_COLUMN, "") if rows else ""


def score_command(docktscore: Path, receptor: Path, pose: Path, cofactors: list[Path], model: str, extra: list[str]):
    """The DockTScore command line that scores one pose with one model."""
    command = [str(docktscore), "-p", str(receptor), "-l", str(pose), "--model", model, "--report-format", "csv"]
    if cofactors:
        command += ["--cofactor", *[str(path) for path in cofactors]]
    return command + extra


def write_scores_csv(path: Path, pdb: str, scores: dict[str, str]) -> None:
    """`<pdb>_scores.csv`: the header and one row, the complex and its score per model column."""
    with open(path, "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["complex", *scores])
        writer.writerow([pdb, *scores.values()])


def docked_poses(out_dir: Path) -> Path | None:
    """The multi-pose file DockThor wrote for the complex, or None."""
    docked = sorted(out_dir.glob("*_docked.mol2"))
    return docked[0] if docked else None


def score_complex(args, pdb: str, out_dir: Path, log_path: Path | None = None) -> dict[str, str]:
    """Score the best-scored pose of a docked complex with every model; {column: score} of the models that ran.

    Needs args.docktscore_path, args.dataset, args.docktscore_models (a list) and args.docktscore_args (a list).
    A complex that cannot be scored at all gives {} after a warning, and the set goes on."""
    docked = docked_poses(out_dir)
    if docked is None:
        print(f"    [warn] {pdb}: no *_docked.mol2 to rescore in {out_dir}", file=sys.stderr)
        return {}
    pose = out_dir / f"{pdb}_top1.mol2"
    if not split_best_pose(docked, pose):
        print(f"    [warn] {pdb}: {docked.name} holds no pose", file=sys.stderr)
        return {}
    try:
        receptor = receptor_pdb(args.dataset, pdb)
    except (NotATopology, OSError) as error:  # a receptor the charged PDB cannot be written from
        print(f"    [warn] {pdb}: no receptor to score against: {error}", file=sys.stderr)
        return {}
    cofactors = cofactor_files(args.dataset, pdb)

    scores = {}
    for model in args.docktscore_models:
        # DockTScore names its output after the ligand, so every model scores its own copy of the pose.
        model_pose = out_dir / f"{pdb}_top1_{model}.mol2"
        shutil.copyfile(pose, model_pose)
        command = score_command(args.docktscore_path, receptor, model_pose, cofactors, model, args.docktscore_args)
        log = log_path or out_dir / f"{pdb}_score.log"
        with open(log, "a") as stream:
            stream.write(f"$ {' '.join(command)}\n")
            stream.flush()
            completed = subprocess.run(command, cwd=out_dir, stdout=stream, stderr=subprocess.STDOUT)
        if completed.returncode != 0:
            print(f"    [warn] {pdb}: docktscore --model {model} failed (exit {completed.returncode})", file=sys.stderr)
            continue
        if score := read_score(out_dir / f"{model_pose.stem}_score.csv"):
            scores[column_of(model)] = score
        else:
            print(f"    [warn] {pdb}: no {SCORE_COLUMN} for model {model}", file=sys.stderr)
    if scores:
        write_scores_csv(out_dir / f"{pdb}{SCORES_CSV}", pdb, scores)
    return scores
