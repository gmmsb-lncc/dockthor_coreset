# dockthor_coreset

Redocking qualification set for the [DockThor Suite](https://github.com/gmmsb-lncc/dockthor-suite): 195 protein–ligand complexes, prepared in DockThor's own formats and ready to dock, with the TOML configuration of every complex and the script that runs the whole set and reports the success rates.

It answers one question about a DockThor build or setting: **how often does DockThor put the ligand back into its crystal pose?** Every complex is redocked into its own receptor, and each GA run is a success when its pose is within a given RMSD of the crystal ligand.

- [1. The dataset](#1-the-dataset)
- [2. Repository layout](#2-repository-layout)
- [3. Quick start](#3-quick-start)
- [4. The qualification script](#4-the-qualification-script)
- [5. Output](#5-output)
- [6. Practical notes](#6-practical-notes)
- [7. References](#7-references)

---

## 1. The dataset

### 1.1 Origin

The complexes are the **core set of PDBbind v2013**, the test set of the CASF-2013 benchmark (Li et al., 2014): 65 protein targets with 3 complexes each, chosen from the PDBbind refined set to cover a wide range of binding affinities with high-quality crystal structures.

Two checks point to v2013 rather than a later release:

- **Size.** The CASF-2013 core set has exactly 195 complexes. The CASF-2016 core set, used with PDBbind v2016 and v2020 (v2020 introduced no new core set), has 285.
- **Dates.** The newest entries here are 2012 depositions (`4g8m`, `4gid`, `4gqq`, `4djv`) and no ID starts with `5`. CASF-2016 includes complexes deposited from 2013 to 2015.

> In this repository the set is also called *vcoreset*. The PDBbind/CASF files were converted to DockThor topologies; the structures are not modified otherwise.

### 1.2 Contents

| | Count | Notes |
|---|---|---|
| Complexes | 195 | Named by PDB ID (`10gs`, `1a30`, ... `4tmn`) |
| Receptors | 195 | DockThor receptor topology (`.in`, the PDBThorBox format) |
| Ligands | 195 | DockThor ligand topology (`.top`, MMFF94S, the MMFFLigand format); 15 to 140 atoms, median 42 |
| Complexes with cofactors | 70 | 239 cofactor topologies in total, up to 23 for one complex (metal ions, waters, coenzymes) |
| Reference poses | 195 | Crystal ligand pose as SDF |
| TOMLs | 195 | One DockThor configuration per complex (`toml/`) |

### 1.3 Files of one complex

For a complex `<pdb>`:

| File | Contents |
|------|----------|
| `receptor/<pdb>.in` | Receptor topology: one line per atom with index, type, coordinates, bonded neighbors, charge, atom/residue name, chain |
| `ligand/<pdb>_ligand.top` | Ligand topology in the **crystal pose**: `$TOTAL_ATOMS` (element, index, MMFF type, charge, x y z, atomic number per atom), `$CON` (bonds), `$SELECTED_TORSIONS` (rotatable bonds) and `$NONBONDED_INTERACTIONS` |
| `cofactor/<pdb>/<pdb>_cofactor-NN.top` | Optional: one topology per cofactor (same format as the ligand), with the matching `.pdb` |
| `config/<pdb>.grid_center.txt` | Center of the docking box, `x y z` in Å |
| `config/<pdb>.grid_size.txt` | Edge of the box, `x y z` in Å; every box is a cube, 11.0 to 26.6 Å (median 16.8) |
| `config/<pdb>.rstep.txt` | Grid step in Å: 0.25 for 186 complexes; for the 9 largest boxes it is raised (up to 0.2987, `3ag9`) to keep 89 intervals per edge, under DockThor's 1,000,000 grid points |
| `reference/<pdb>_ligand.sdf` | Crystal ligand pose (converted with Open Babel), for external RMSD checks |
| `toml/<pdb>.toml` | DockThor configuration of the complex (§1.4) |

A few extra files are kept for reference: `ligand/<pdb>_ligand.top_diedral.inf` (rotatable dihedrals written by the ligand preparation, for `10gs`, `1n1m`, `4dew`, `4tmn`) and `reference/10gs_ligand.top` / `10gs_ligand_with_H.pdb`. DockThor does not read them.

**The reference of the success rate is the ligand topology itself.** DockThor measures the RMSD of each docked pose against the coordinates of the input `.top`, which are the crystal pose. The SDFs in `reference/` hold the same heavy-atom positions (checked for all 195 complexes) and are there for tools outside DockThor (`obrms`, RDKit).

### 1.4 The TOMLs

Each `toml/<pdb>.toml` holds the inputs and the box of one complex, nothing about the docking itself:

```toml
## DockThor qualification set: 1f8b

output-dir="vcoreset_run"
ligand="ligand/1f8b_ligand.top"
receptor="receptor/1f8b.in"
cofactors-path="cofactor/1f8b"          # only for the 70 complexes with cofactors
grid-size=[14.022,14.022,14.022]
grid-center=[26.464,17.2705,62.4755]
grid-rstep=0.25
```

- **Paths are relative to the repository root**, so DockThor must run from there; the script does that (§4.3). To run one TOML by hand: `cd dockthor_coreset; dockthor --config toml/1f8b.toml`.
- **The docking setup is not in the TOMLs.** The evolution type, number of runs, success-rate cutoff and the rest are given when the set is run (`--dockthor-args`), so the same TOMLs serve every setup.
- `output-dir` is `vcoreset_run`; the script's `-o` overrides it.

---

## 2. Repository layout

```
dockthor_coreset/
├── README.md
├── scripts/
│   └── qualification_set.py   generate and run the TOMLs, collect the success rates
├── toml/                      195 TOMLs, one per complex
├── receptor/                  <pdb>.in
├── ligand/                    <pdb>_ligand.top
├── cofactor/                  <pdb>/<pdb>_cofactor-NN.{top,pdb}   (70 complexes)
├── config/                    <pdb>.grid_center.txt, <pdb>.grid_size.txt, <pdb>.rstep.txt
└── reference/                 <pdb>_ligand.sdf
```

Results are written inside the repository, one directory per setup (`steadystate/`, `generational_ocl/`, ...); they are untracked files, not part of the dataset.

---

## 3. Quick start

Requirements: Python 3.11+ (standard library only), git, and a DockThor build from [dockthor-suite](https://github.com/gmmsb-lncc/dockthor-suite) (its `README.md` explains the build; `build/bin/dockthor` for the CPU, `build_ocl/bin/dockthor` for OpenCL). The repository is private: clone it with SSH or your git credentials.

```fish
git clone git@github.com:gmmsb-lncc/dockthor_coreset.git
cd dockthor_coreset
set DOCKTHOR ~/dockthor-suite        # your dockthor-suite checkout (bash: DOCKTHOR=~/dockthor-suite)
```

**SteadyState on the CPU, MPI.** 24 GA runs per complex over 24 MPI processes (one run each), success rate up to 3 Å:

```fish
python3 scripts/qualification_set.py --run-dir toml -o steadystate \
    --run-path $DOCKTHOR/build/bin/dockthor --mpi --np 24 --skip-done \
    --dockthor-args='--ga-evolution-type STEADYSTATE --ga-nrun 24 --success-rate-RMSD 3'
```

**GENERATIONAL_OCL on the GPU.** OpenCL build, one process (the GPU packs the 24 runs, so no `--mpi`):

```fish
python3 scripts/qualification_set.py --run-dir toml -o generational_ocl \
    --run-path $DOCKTHOR/build_ocl/bin/dockthor --skip-done \
    --dockthor-args='--ga-evolution-type GENERATIONAL_OCL --ga-nrun 24 --success-rate-RMSD 3'
```

Each prints one line per complex and ends with the summary path:

```
Dataset root (from the TOMLs): /home/user/dockthor_coreset
[1/195] 10gs: mpirun -n 24 .../dockthor --config .../toml/10gs.toml --output-dir steadystate --ga-evolution-type STEADYSTATE --ga-nrun 24 --success-rate-RMSD 3
    ok in <t>s, log: .../steadystate/10gs/10gs_run.log
...
Ran 195 system(s) in <t>s, 0 failed.
Summary: /home/user/dockthor_coreset/steadystate/qualification_summary.csv
```

`steadystate/qualification_summary.csv` then holds the execution time and both success rates of every complex (§5.3). If the run stops, the same command resumes it: `--skip-done` skips the finished complexes.

---

## 4. The qualification script

`scripts/qualification_set.py` does three things, in any combination:

1. **Generate** a TOML per complex from `config/`, `ligand/`, `receptor/` and `cofactor/` (`--generate`).
2. **Run** DockThor on a folder of TOMLs, saving each complex's console output (`--run-path`).
3. **Summarize** the execution time and success rates of every complex run (automatic after 2).

The dataset can be this repository, any directory with the same layout, or a git repository the script clones (`--git-url`).

### 4.1 Workflows

**A. Run the committed TOMLs** (the usual case, §3). `--run-dir toml` takes every `toml/<pdb>.toml`:

```fish
python3 scripts/qualification_set.py --run-dir toml -o steadystate --run-path $DOCKTHOR/build/bin/dockthor \
    --mpi --np 24 --dockthor-args='--ga-nrun 24 --success-rate-RMSD 3'
```

**B. A subset.** `--systems` restricts any workflow to some complexes, for a quick check of a build:

```fish
python3 scripts/qualification_set.py --run-dir toml --systems 1bcu 10gs 1f8b -o check \
    --run-path $DOCKTHOR/build/bin/dockthor --mpi --np 4 --dockthor-args='--ga-nrun 4 --success-rate-RMSD 2'
```

**C. Generate TOMLs with the setup inside them.** `--generate` writes new TOMLs, `--set` adds options to all of them, and `--output-toml-path` keeps them apart from `toml/`. The setup is then recorded with the TOMLs, and `--dockthor-args` is not needed:

```fish
python3 scripts/qualification_set.py --generate -o generational_ocl --output-toml-path toml_generational_ocl \
    --set 'ga-evolution-type="GENERATIONAL_OCL"' --set ga-nrun=24 --set success-rate-RMSD=3.0 \
    --run-path $DOCKTHOR/build_ocl/bin/dockthor --skip-done
```

Without `--output-toml-path` the TOMLs are written to the dataset root; existing ones are kept unless `--overwrite` is given.

**D. Straight from GitHub.** `--git-url` clones the set (into `./dockthor_coreset` by default; an existing clone there is reused) and runs any of the above on the clone; `--git-ref` pins a branch, tag or commit:

```fish
python3 dockthor_coreset/scripts/qualification_set.py --git-url git@github.com:gmmsb-lncc/dockthor_coreset.git \
    --git-ref main --run-dir toml -o steadystate --run-path $DOCKTHOR/build/bin/dockthor --mpi --np 24 \
    --dockthor-args='--ga-nrun 24 --success-rate-RMSD 3'
```

The script prints the checked-out commit (`Dataset commit: <sha>`), which ties the results to an exact version of the set.

**E. Another dataset.** Any directory with `config/`, `ligand/`, `receptor/` (and optionally `cofactor/`) laid out as in §1.3 works; pass it as the first argument:

```fish
python3 scripts/qualification_set.py ~/my_set --generate -o results --set ga-nrun=10 \
    --run-path $DOCKTHOR/build/bin/dockthor
```

### 4.2 Options

Values that start with `-` must use the `=` form: `--dockthor-args='--seed 7'` works, `--dockthor-args '--seed 7'` does not.

**Dataset**

| Option | Default | Description |
|--------|---------|-------------|
| `dataset` (positional) | this repository | Root of the set. Omit it to use the repository holding the script, with `--git-url`, or with `--run-dir` (the root is then inferred from the TOMLs, §4.3) |
| `--systems PDB ...` | all | Only these complexes. Without it: every TOML of `--run-dir`, or every complex of `config/` |
| `--git-url URL` | — | Clone the set from a git repository (HTTPS, SSH or a local path) instead of giving `dataset` |
| `--git-ref REF` | default branch | Branch, tag or commit to check out. Branches and tags are cloned shallow (`--depth 1`); a commit needs a full clone, done automatically |
| `--git-subdir PATH` | repository root | Set root inside the cloned repository, when the set is in a subfolder |
| `--clone-dir DIR` | `./<repo name>` | Where to clone. An existing clone there is reused as is; a non-empty folder that is not a clone stops the script |
| `--git-update` | off | Fetch `--git-ref` (or the remote default branch) into an existing clone and check it out before running |

**Generation**

| Option | Default | Description |
|--------|---------|-------------|
| `--generate` | off | Write `<pdb>.toml` for every selected complex. Requires `-o` and a dataset with `config/` |
| `-o, --output-dir DIR` | TOML's `output-dir` | With `--generate`, the `output-dir` written into the TOMLs. With `--run-path`, also passed to DockThor as `--output-dir`, overriding the TOMLs, so the same TOMLs can be run into one directory per setup. A relative path is relative to the dataset root |
| `--output-toml-path DIR` | dataset root | Where `--generate` writes the TOMLs and where `--run-path` reads them (when `--run-dir` is not given). Relative to the current directory. Input paths inside the TOMLs stay relative to the dataset root |
| `--set KEY=VALUE` | — | Any DockThor option, written into every generated TOML. `VALUE` is a TOML literal: numbers and booleans as is (`--set ga-nrun=24`, `--set symmetric-RMSD=true`), strings quoted (`--set 'ga-evolution-type="STEADYSTATE"'`). Repeatable |
| `--no-cofactors` | off | Leave `cofactors-path` out of the generated TOMLs, to dock without cofactors |
| `--overwrite` | off | Replace existing TOMLs (otherwise they are kept and reported as skipped) |

**Execution**

| Option | Default | Description |
|--------|---------|-------------|
| `--run-path DOCKTHOR` | — | The `dockthor` binary. Runs every selected TOML, one complex after the other |
| `--run-dir DIR` | — | Run the TOMLs already in `DIR` instead of generated ones. Alone, `DIR` is relative to the current directory; with `dataset` or `--git-url`, relative to the dataset root |
| `--mpi` | off | Launch DockThor with `mpirun` |
| `--np N` | 4 | MPI processes per complex. DockThor spreads the `ga-nrun` runs over them, so `--np` equal to `ga-nrun` runs them all at once |
| `--mpirun CMD` | `mpirun` | MPI launcher, e.g. `srun` on a SLURM cluster |
| `--mpi-args=ARGS` | — | Extra launcher arguments, e.g. `--mpi-args='--oversubscribe --bind-to none'` |
| `--dockthor-args=ARGS` | — | Extra DockThor command-line options for every complex; they take precedence over the TOMLs |
| `--summary-csv FILE` | §5.3 | Path of the summary CSV |
| `--skip-done` | off | Skip complexes that already have a `*_docked.csv` in `<output-dir>/<pdb>/`: rerun the same command to resume |
| `--stop-on-error` | off | Stop at the first complex whose DockThor call fails (default: report it and go on) |
| `--dry-run` | off | Print the DockThor commands without running them |
| `-v, --verbose` | off | Also print DockThor's output to the terminal (it is always saved to the run log) |

The exit status is non-zero when any complex fails or a TOML cannot be generated.

**Useful DockThor options for `--dockthor-args` or `--set`** (the full list is in the DockThor manual, `doc/MANUAL_dmrts.md` §6, or `dockthor --help`):

| Option | Meaning |
|--------|---------|
| `--ga-evolution-type` | `STEADYSTATE` (CPU, default) or `GENERATIONAL_OCL` (GPU, OpenCL build) |
| `--ga-nrun N` | Independent GA runs per complex: the denominator of the success rate. 24 runs give 4.2% steps |
| `--success-rate-RMSD c` | Turn the success rate on and report it at 0, 0.5, ... up to `c` Å. **Needed for the summary rates** |
| `--symmetric-RMSD` | Symmetry-aware RMSD for clustering and the success rate. More accurate for symmetric ligands, much slower (§6) |
| `--seed N` | RNG seed, for repeating a qualification exactly |
| `--ga-evaluations N`, `--ga-init-evaluations N` | Length of the evolution (lower for quick checks) |
| `--ocl-device-type`, `--ocl-platform-name` | GPU selection for `GENERATIONAL_OCL`; `dockthor --ocl-display-devices` lists the devices |

### 4.3 How paths are resolved

The input paths in the TOMLs (`ligand/...`, `receptor/...`, `cofactor/...`) are relative to the **dataset root**, and the script always runs DockThor with the dataset root as working directory, passing each TOML by its absolute path. The dataset root is:

| Given | Dataset root | `--run-dir DIR` relative to |
|-------|--------------|-----------------------------|
| nothing (only `--run-path`, `--generate`, ...) | the repository holding the script | the dataset root |
| `dataset` | `dataset` | the dataset root |
| `--git-url` | the clone (or `--git-subdir` in it) | the dataset root |
| `--run-dir DIR` only | inferred: `DIR` or the first parent of it where the `ligand` and `receptor` of a TOML exist | the current directory |

With `--run-dir` alone the script prints the root it inferred (`Dataset root (from the TOMLs): ...`) and stops if no such folder exists. Run the examples of this README from the repository root, or give `--run-dir` as an absolute path.

The output directory (`-o` or the TOML's `output-dir`) is relative to the dataset root too: `-o steadystate` writes to `dockthor_coreset/steadystate/`.

---

## 5. Output

### 5.1 Directory tree

```
<output-dir>/
├── qualification_summary.csv                    all complexes run (§5.3)
├── 10gs/
│   ├── 10gs_run.log                             console output of the DockThor call
│   ├── 10gs_ligand_docked.csv / .mol2           cluster leaders of all runs, best energy first
│   ├── 10gs_ligand_docked_success_rate.csv      success rate per cutoff (with --success-rate-RMSD)
│   ├── 10gs_ligand_docked_success_t_energy.csv / .mol2   best-energy pose of every run
│   └── 10gs_ligand_docked_success_RMSD.csv / .mol2       pose closest to the crystal of every run
├── 1a30/
...
```

Each complex directory is named after its receptor file (the PDB ID). The DockThor files are described in the DockThor manual (`doc/MANUAL_dmrts.md` §10).

### 5.2 Run log

`<output-dir>/<pdb>/<pdb>_run.log` is the console output of the DockThor call: the command line, then stdout and stderr, without the terminal colors. It holds the grid, the evolution progress and the success-rate table:

```
[DOCKTHOR] Success rate (%) over 24 runs, RMSD <= cutoff.
           |------->     RMSD     energy       RMSD
           |------->     0.00       0.00       0.00
           |------->     0.50      25.00      25.00
           |------->     1.00      33.33      41.67
           ...
```

### 5.3 Summary CSV

After every complex the script rewrites `qualification_summary.csv`, so an interrupted set keeps the rows already done. It goes to the common output directory (`<output-dir>/qualification_summary.csv`), or next to the TOMLs when the complexes write to different directories; `--summary-csv FILE` sets it explicitly.

```
complex,execution_time,RMSD_0,RMSD_0.5,RMSD_1,RMSD_1.5,RMSD_2,RMSD_2.5,RMSD_3,Energy_0,Energy_0.5,Energy_1,Energy_1.5,Energy_2,Energy_2.5,Energy_3
10gs,248.70,0.00,12.50,37.50,41.67,41.67,41.67,41.67,0.00,8.33,25.00,33.33,33.33,37.50,41.67
1a30,...
```

| Column | Unit | Description |
|--------|------|-------------|
| `complex` | — | PDB ID |
| `execution_time` | s | Wall time of the DockThor call (with `mpirun`), measured by the script |
| `RMSD_<c>` | % | Success rate **by RMSD** at cutoff `c` Å |
| `Energy_<c>` | % | Success rate **by energy** at cutoff `c` Å |

The rates are filled only when DockThor runs with `--success-rate-RMSD`; otherwise the script warns and the CSV has the times only. A failed complex keeps its row with empty rates; complexes skipped by `--skip-done` are not listed (they were not run in that call). When complexes use different cutoffs, the columns are the union and missing rates are empty.

### 5.4 What the success rates mean

DockThor computes the rates at the end of each docking; the full definition is in the DockThor manual, `doc/MANUAL_dmrts.md` §10.5. In short, with *N* = `ga-nrun` runs:

1. Each run is clustered on its own, and its leaders are sorted by energy.
2. **By energy:** each run contributes its lowest-energy leader, the pose DockThor would answer with.
3. **By RMSD:** each run contributes its leader closest to the crystal pose, which shows whether the run found the native pose at all.
4. A run succeeds at cutoff `c` when that pose's heavy-atom RMSD to the crystal ligand is `<= c` Å; the rate is `successful runs / N × 100`.

The rate by RMSD is never below the rate by energy; the gap is the share of runs that found the native pose but scored another one better. The usual redocking criterion is **2.0 Å**. To pool the set, average a column over the complexes (each complex weighs the same), for example:

```fish
python3 -c "import csv,sys; r=list(csv.DictReader(open(sys.argv[1]))); c=sys.argv[2]; v=[float(x[c]) for x in r if x[c]]; print(f'{c}: {sum(v)/len(v):.2f}% over {len(v)} complexes')" steadystate/qualification_summary.csv Energy_2
```

---

## 6. Practical notes

- **One output directory per setup.** Give every setup its own `-o` (`steadystate`, `generational_ocl`, `steadystate_sym`, ...). With `--grid-write`, DockThor reuses a grid file already in the output directory after checking only its size, receptor atom count, center and atom types, so a grid left from another build or setting can be loaded silently. Use a fresh directory, or do not pass `--grid-write`.
- **MPI processes.** With `--np` equal to `ga-nrun` all runs evolve at once; with fewer processes each one runs several in turn (`--np 12` with 24 runs takes about twice as long per complex). `GENERATIONAL_OCL` needs a single process: no `--mpi`.
- **Open MPI accelerator conflict.** If `mpirun` fails with a rocm/cuda accelerator error, export `OMPI_MCA_accelerator='^cuda,rocm'` (fish: `set -x OMPI_MCA_accelerator '^cuda,rocm'`) before running.
- **Symmetric RMSD is slow.** `--symmetric-RMSD` gives lower, more correct RMSDs for ligands with symmetric groups (on `10gs`, 0.91 → 0.38 Å for one pose), which raises the rates at tight cutoffs; at 2 Å it rarely changes them. It made the clustering of `10gs` take about 6.5 minutes instead of seconds, so expect a much longer set.
- **Resuming.** `--skip-done` skips complexes with results, so a killed set resumes with the same command; the summary of the new call lists only the complexes it ran. Delete a complex's directory to run it again.
- **Clusters.** Use `--mpirun srun` (and `--mpi-args` for its options) inside a SLURM allocation; the run logs and the summary are written as the set progresses.
- **Checking the RMSDs outside DockThor.** `reference/<pdb>_ligand.sdf` with Open Babel's symmetry-aware `obrms` reproduces DockThor's `--symmetric-RMSD` values: `obrms -f reference/10gs_ligand.sdf steadystate/10gs/10gs_ligand_docked.mol2`.

---

## 7. References

- **CASF-2013 core set:** Li, Y.; Liu, Z.; Li, J.; Han, L.; Liu, J.; Zhao, Z.; Wang, R. Comparative Assessment of Scoring Functions on an Updated Benchmark: 1. Compilation of the Test Set. *J. Chem. Inf. Model.* **2014**, 54, 1700–1716. DOI [10.1021/ci500080q](https://doi.org/10.1021/ci500080q).
- **CASF-2013 evaluation:** Li, Y.; Han, L.; Liu, Z.; Wang, R. Comparative Assessment of Scoring Functions on an Updated Benchmark: 2. Evaluation Methods and General Results. *J. Chem. Inf. Model.* **2014**, 54, 1717–1736. DOI [10.1021/ci500081m](https://doi.org/10.1021/ci500081m).
- **PDBbind:** Wang, R.; Fang, X.; Lu, Y.; Wang, S. The PDBbind Database: Collection of Binding Affinities for Protein–Ligand Complexes with Known Three-Dimensional Structures. *J. Med. Chem.* **2004**, 47, 2977–2980. DOI [10.1021/jm030580l](https://doi.org/10.1021/jm030580l).
- **DockThor:** see the references of the [DockThor Suite](https://github.com/gmmsb-lncc/dockthor-suite) (`doc/MANUAL_dmrts.md` §16).

The PDBbind structures are distributed under the PDBbind terms of use; cite PDBbind and CASF when using this set.
