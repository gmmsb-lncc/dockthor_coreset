#!/usr/bin/env bash
# Dock one complex of the set with one setup inside the private DockThor release image, as a CI job does.
#
#   scripts/ci/dock.sh <setup> <pdb> <image>
#   scripts/ci/dock.sh steadystate_symmetric 1bcu ghcr.io/gmmsb-lncc/dockthor:2.1.0
#
# Run from the repository root. Results go to ci_out/<setup>/<pdb>/; the job keeps only
# <pdb>_ligand_docked_success_rate.csv (scripts/qualification_set.py --collect merges them).
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

# The setups of the release qualification (README §7.1). Both use the symmetry-aware RMSD to the crystal ligand.
common_args="--ga-nrun 24 --success-rate --success-rate-by-run 3 --symmetric-RMSD"
case "$setup" in
    steadystate_symmetric)
        binary=dockthor-mpi
        # 24 runs over 2 MPI processes. --oversubscribe: Open MPI counts physical cores as slots, and a runner's
        # container may show fewer than its vCPUs.
        launch=(--mpi --np 2 --mpi-args=--oversubscribe)
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

if [ ! -f "toml/$pdb.toml" ]; then
    echo "no toml/$pdb.toml" >&2
    exit 2
fi

# Own user: outputs belong to the runner user, and mpirun refuses root. HOME=/tmp: writable caches (Open MPI, PoCL).
docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp -v "$PWD:/data" -w /data "$image" \
    python3 scripts/qualification_set.py --config-dir toml --systems "$pdb" -o "ci_out/$setup" \
        --dockthor-path "/usr/local/bin/$binary" "${launch[@]}" --dockthor-args="$dockthor_args"
