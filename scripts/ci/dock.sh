#!/usr/bin/env bash
# Dock one complex of the set with one setup inside the private DockThor release image, as a CI job does.
#
#   scripts/ci/dock.sh <setup> <pdb> <image>
#   scripts/ci/dock.sh steadystate_symmetric 1bcu ghcr.io/gmmsb-lncc/dockthor:2.1.0
#
# SET_DIR picks the qualification set (default v2013-coreset); the repository can hold several.
#
# Run from the repository root. Results go to ci_out/<setup>/<pdb>/ at the repository root -- not inside the set,
# which is why -o is given as an absolute container path (-o is otherwise relative to the set root); the job keeps only
# <pdb>_ligand_docked_success_rate.csv and <pdb>_scores.csv, the row of the rescoring
# (scripts/summarize_results.py merges both).
# DockThor's console output stays in ci_out/<setup>/<pdb>/<pdb>_run.log on the runner and is never printed:
# this repository and its Actions logs are public.
set -euo pipefail

if [ $# -ne 3 ]; then
    echo "usage: $0 <setup> <pdb> <image>" >&2
    exit 2
fi
setup="$1"
pdb="$2"
image="$3"
set_dir="${SET_DIR:-v2013-coreset}"

# The setups of the release qualification (README §7.1). Both use the symmetry-aware RMSD to the crystal ligand.
common_args="--ga-nrun 24 --success-rate --success-rate-by-run 3 --symmetric-RMSD"
case "$setup" in
    steadystate_symmetric)
        binary=dockthor-mpi
        # 24 runs over 4 MPI processes, one per vCPU of a public repository's runner. --oversubscribe: Open MPI counts
        # physical cores as slots, and a runner's container may show fewer than its vCPUs.
        launch=(--mpi --np 4 --mpi-args=--oversubscribe)
        dockthor_args="--ga-evolution-type STEADYSTATE $common_args"
        ;;
    generational_ocl_symmetric)
        binary=dockthor-ocl
        # One process: the OpenCL device (PoCL on the runner's CPU) packs the 24 runs. No MPI.
        launch=()
        dockthor_args="--ga-evolution-type GENERATIONAL_OCL $common_args"
        ;;
    *)
        echo "unknown setup: $setup (steadystate_symmetric, generational_ocl_symmetric)" >&2
        exit 2
        ;;
esac

if [ ! -f "$set_dir/$pdb/ligands/${pdb}_ligand.top" ] || [ ! -f "$set_dir/$pdb/receptors/$pdb.in" ]; then
    echo "no $set_dir/$pdb/ligands/${pdb}_ligand.top or $set_dir/$pdb/receptors/$pdb.in" >&2
    exit 2
fi

# Own user: outputs belong to the runner user, and mpirun refuses root. HOME=/tmp: writable caches (Open MPI, PoCL).
in_image=(docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp -v "$PWD:/data" -w /data "$image")

# The TOMLs are generated from the ligand, not committed (README §1.5): one complex takes milliseconds. They hold
# absolute paths, so they are generated in the container, where the repository is /data.
"${in_image[@]}" python3 scripts/grid_toml_generator.py "$set_dir" --complex "$pdb" --overwrite

# DockTScore rescores the best-scored pose in the same run, which fills the affinity statistics and their badges
# (README 4.7). A release image without docktscore simply docks: the docking results are the same either way.
rescore=()
if "${in_image[@]}" test -x /usr/local/bin/docktscore; then
    rescore=(--docktscore-path /usr/local/bin/docktscore)
else
    echo "note: no docktscore in $image, so no rescoring (the affinity badges keep their last values)" >&2
fi

"${in_image[@]}" python3 scripts/qualification_set.py --config-dir "$set_dir/config" --complex "$pdb" -o "/data/ci_out/$setup" \
    --dockthor-path "/usr/local/bin/$binary" "${launch[@]}" --dockthor-args="$dockthor_args" "${rescore[@]}"
