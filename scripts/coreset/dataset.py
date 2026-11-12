"""The qualification set on disk: its systems, the boxes of config/, the cofactors, and cloning it from git.

Layout of a set root:

    config/<pdb>.grid_center.txt   three floats (x y z)
    config/<pdb>.grid_size.txt     three floats (x y z), in Angstroms
    config/<pdb>.rstep.txt         one float (grid step)
    ligand/<pdb>_ligand.top        ligand topology
    receptor/<pdb>.in              receptor topology
    cofactor/<pdb>/                optional, cofactor topologies
"""

import re
import shlex
import subprocess
import sys
import tomllib
from pathlib import Path


def read_floats(path: Path, count: int) -> list[float]:
    values = [float(v) for v in path.read_text().split()]
    if len(values) != count:
        raise ValueError(f"{path}: expected {count} value(s), found {len(values)}")
    return values


def git(*command: str, cwd: Path | None = None) -> subprocess.CompletedProcess:
    print("    $ git " + shlex.join(command), flush=True)
    return subprocess.run(["git", "-c", "advice.detachedHead=false", *command], cwd=cwd)


def repo_name(url: str) -> str:
    """https://github.com/org/set.git, git@github.com:org/set.git or /path/set -> set"""
    return re.split(r"[/:]", url.rstrip("/"))[-1].removesuffix(".git")


def clone_dataset(args) -> Path:
    """Clone (or reuse) the repository of --git-url and return the set root inside it."""
    clone_dir = (args.clone_dir or Path(repo_name(args.git_url))).expanduser().resolve()
    if (clone_dir / ".git").is_dir():
        print(f"Using the existing clone {clone_dir}")
        if args.git_update:
            ref = args.git_ref or "HEAD"
            if (
                git("fetch", "--depth", "1", "origin", ref, cwd=clone_dir).returncode != 0
                or git("checkout", "--detach", "FETCH_HEAD", cwd=clone_dir).returncode != 0
            ):
                sys.exit(f"[error] could not update {clone_dir} to {ref}")
    elif clone_dir.exists() and any(clone_dir.iterdir()):
        sys.exit(f"[error] {clone_dir} exists and is not a git clone; choose another --clone-dir")
    else:
        print(f"Cloning {args.git_url} into {clone_dir}")
        branch = ["--branch", args.git_ref] if args.git_ref else []
        if git("clone", "--depth", "1", *branch, args.git_url, str(clone_dir)).returncode != 0:
            if not args.git_ref:
                sys.exit(f"[error] could not clone {args.git_url}")
            # --branch takes only branches and tags: a commit needs a full clone and a checkout.
            if clone_dir.exists() and any(clone_dir.iterdir()):
                sys.exit(f"[error] could not clone {args.git_url}; {clone_dir} was left behind, remove it")
            if (
                git("clone", args.git_url, str(clone_dir)).returncode != 0
                or git("checkout", "--detach", args.git_ref, cwd=clone_dir).returncode != 0
            ):
                sys.exit(f"[error] could not clone {args.git_url} at {args.git_ref}")
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=clone_dir, capture_output=True, text=True).stdout.strip()
    print(f"Dataset commit: {commit}")
    return clone_dir / args.git_subdir if args.git_subdir else clone_dir


def infer_dataset(toml_dir: Path) -> Path | None:
    """The set root of a TOML folder: the first of toml_dir and its parents where the ligand and
    receptor of a TOML exist (their paths are relative to the set root)."""
    for toml_path in sorted(toml_dir.glob("*.toml")):
        try:
            config = tomllib.loads(toml_path.read_text())
        except (OSError, tomllib.TOMLDecodeError):
            continue
        if "ligand" not in config or "receptor" not in config:
            continue
        for root in (toml_dir, *toml_dir.parents):
            if (root / config["ligand"]).exists() and (root / config["receptor"]).exists():
                return root
        return None
    return None


def discover_systems(root: Path) -> list[str]:
    return sorted(p.name.removesuffix(".grid_center.txt") for p in (root / "config").glob("*.grid_center.txt"))


def cofactors_dir(root: Path, cofactor_path: Path, pdb: str) -> Path | None:
    """cofactor_path/<pdb> (relative to root, as written in a TOML) when it holds cofactor topologies, else None."""
    cofactors = cofactor_path / pdb
    return cofactors if any((root / cofactors).glob("*.top")) else None


def build_config(root: Path, pdb: str, output_dir: Path, use_cofactors: bool) -> dict:
    """Collect the DockThor options of one system from config/*.txt. Input paths are relative to the set root."""
    config_dir = root / "config"
    ligand = Path("ligand") / f"{pdb}_ligand.top"
    receptor = Path("receptor") / f"{pdb}.in"
    for path in (ligand, receptor):
        if not (root / path).is_file():
            raise FileNotFoundError(f"missing {root / path}")

    config = {
        "output-dir": str(output_dir),
        "ligand": str(ligand),
        "receptor": str(receptor),
    }
    if use_cofactors and (cofactors := cofactors_dir(root, Path("cofactor"), pdb)) is not None:
        config["cofactors-path"] = str(cofactors)
    config["grid-size"] = read_floats(config_dir / f"{pdb}.grid_size.txt", 3)
    config["grid-center"] = read_floats(config_dir / f"{pdb}.grid_center.txt", 3)
    config["grid-rstep"] = read_floats(config_dir / f"{pdb}.rstep.txt", 1)[0]
    return config
