#!/usr/bin/env bash
# Print, as a JSON list, the complexes of $SYSTEMS (a JSON list of PDB IDs) that have no success-rate report
# artifact sr-<setup>-<pdb> in this workflow run yet: the ones whose dock job failed (qualification.yml docks them
# again before collecting).
#
#   SYSTEMS='["1bcu","2vw5"]' GH_TOKEN=... scripts/ci/missing_reports.sh <setup>
set -euo pipefail

if [ $# -ne 1 ]; then
    echo "usage: $0 <setup>" >&2
    exit 2
fi
setup="$1"

names=$(gh api --paginate "repos/$GITHUB_REPOSITORY/actions/runs/$GITHUB_RUN_ID/artifacts?per_page=100" \
    --jq '.artifacts[] | select(.expired | not) | .name')

jq -c --arg setup "$setup" --arg names "$names" \
    '($names | split("\n")) as $have | map(select((("sr-" + $setup + "-" + .) as $n | $have | index($n)) == null))' \
    <<<"$SYSTEMS"
