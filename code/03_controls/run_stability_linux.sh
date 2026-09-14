#!/usr/bin/env bash
# =============================================================================
# Run the stability replicates on any Linux/macOS machine (no scheduler needed).
#
#   ./run_stability_linux.sh [parallel] [replicates] [genes] [cells]
#   defaults: 8 10 800 800
#
# Each job is single-threaded and peaks at ~2 GB RAM (800 genes). Choose
# `parallel` so that parallel x 2 GB stays below available memory.
# Jobs whose output file already exists are skipped, so the script is safe to
# re-run after an interruption.
# =============================================================================
set -euo pipefail

PAR=${1:-8}
REPS=${2:-10}
GENES=${3:-800}
CELLS=${4:-800}

mkdir -p logs work/stability

echo "parallel=${PAR} replicates=${REPS} genes=${GENES} cells=${CELLS}"
echo "started $(date -Is)"

{
  for rep in $(seq 1 "${REPS}"); do
    echo "MSvC_Control ${rep}"
  done
  for rep in $(seq 1 "${REPS}"); do
    echo "MSvC_MS ${rep}"
  done
} | xargs -P "${PAR}" -n 2 bash -c '
  tag="$0"; rep="$1"
  log="logs/stability_${tag}_rep${rep}.log"
  Rscript code/03_controls/28_stability_replicates.R "${tag}" "${rep}" '"${CELLS}"' '"${GENES}"' \
      > "${log}" 2>&1
  echo "finished ${tag} rep ${rep} $(date -Is)"
'

echo "all jobs finished $(date -Is)"
echo
echo "Now aggregate:"
echo "  Rscript code/03_controls/29_stability_analysis.R"
