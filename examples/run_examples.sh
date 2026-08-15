#!/usr/bin/env bash
#
# Run every way of reaching DevKidCC, over every input format, and report.
#
#   bash run_examples.sh <work_dir> [dataset ...]
#
# Expects to run **inside the container**, with <work_dir> holding, for each
# dataset name, a <name>.h5ad and a <name>.rds built from the same cells.
# From the host:
#
#   singularity exec --bind $PWD:/data dkcc.sif \
#       bash /data/examples/run_examples.sh /data organoid_howden fetal_menon
#
# The matrix, per dataset:
#
#   .h5ad  CLI          /opt/dkcc  -> routes to /opt/run_dkcc.py
#   .h5ad  interactive  examples/classify_h5ad.py
#   .rds   CLI          /opt/dkcc  -> routes to /opt/run_dkcc.R
#   .rds   interactive  examples/classify_seurat.R
#
# Then the two interactive outputs are compared cell by cell, since both formats
# hold the same cells. Four runs and one comparison per dataset.

set -uo pipefail

WORK="${1:?usage: run_examples.sh <work_dir> [dataset ...]}"
shift
DATASETS=("$@")
[ ${#DATASETS[@]} -eq 0 ] && DATASETS=(organoid_howden fetal_menon)

EX="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS="$WORK/results"
mkdir -p "$RESULTS"

pass=0; fail=0
declare -a SUMMARY

run_step() {
    local label="$1"; shift
    echo
    echo "--------------------------------------------------------------"
    echo ">>> $label"
    echo "--------------------------------------------------------------"
    local start; start=$(date +%s)
    if "$@" 2>&1; then
        local elapsed=$(( $(date +%s) - start ))
        echo "[PASS] $label (${elapsed}s)"
        pass=$((pass+1)); SUMMARY+=("PASS  ${elapsed}s  $label")
    else
        local elapsed=$(( $(date +%s) - start ))
        echo "[FAIL] $label (${elapsed}s)"
        fail=$((fail+1)); SUMMARY+=("FAIL  ${elapsed}s  $label")
    fi
}

echo "=============================================================="
echo "DevKidCC example run"
echo "  work dir : $WORK"
echo "  examples : $EX"
echo "  datasets : ${DATASETS[*]}"
echo "=============================================================="

for ds in "${DATASETS[@]}"; do
    h5ad="$WORK/$ds.h5ad"
    rds="$WORK/$ds.rds"

    for f in "$h5ad" "$rds"; do
        [ -f "$f" ] || { echo "MISSING INPUT: $f"; fail=$((fail+1)); }
    done

    # 1. h5ad through the CLI dispatcher; must route to Python.
    run_step "$ds  .h5ad  CLI (/opt/dkcc)" \
        /opt/dkcc --input "$h5ad" --output "$RESULTS/${ds}_cli.h5ad"

    # 2. h5ad through the interactive Python example.
    run_step "$ds  .h5ad  interactive (classify_h5ad.py)" \
        python "$EX/classify_h5ad.py" \
            --input "$h5ad" \
            --output "$RESULTS/${ds}_py.h5ad" \
            --summary "$RESULTS/${ds}_py_summary.csv"

    # 3. rds through the CLI dispatcher; must route to R.
    run_step "$ds  .rds   CLI (/opt/dkcc)" \
        /opt/dkcc --input "$rds" --output "$RESULTS/${ds}_cli.rds"

    # 4. rds through the interactive R example.
    run_step "$ds  .rds   interactive (classify_seurat.R)" \
        Rscript "$EX/classify_seurat.R" \
            --input "$rds" \
            --output "$RESULTS/${ds}_r.rds" \
            --summary "$RESULTS/${ds}_r_summary.csv"

    # 5. Same cells, both routes: do they agree?
    if [ -f "$RESULTS/${ds}_r.rds" ] && [ -f "$RESULTS/${ds}_py.h5ad" ]; then
        Rscript "$EX/export_labels.R" \
            "$RESULTS/${ds}_r.rds" "$RESULTS/${ds}_r_labels.csv" >/dev/null 2>&1
        run_step "$ds  route agreement (Python h5ad vs R rds)" \
            python "$EX/compare_routes.py" \
                --h5ad "$RESULTS/${ds}_py.h5ad" \
                --labels "$RESULTS/${ds}_r_labels.csv" \
                --out "$RESULTS/${ds}_agreement.csv"
    fi
done

echo
echo "=============================================================="
echo "SUMMARY"
echo "=============================================================="
for line in "${SUMMARY[@]}"; do echo "  $line"; done
echo
echo "  $pass passed, $fail failed"
echo "  outputs in $RESULTS"

[ "$fail" -eq 0 ] || exit 1
