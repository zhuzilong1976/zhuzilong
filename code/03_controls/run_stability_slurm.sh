#!/bin/bash
# =============================================================================
# Run the stability replicates on a SLURM cluster.
#
# 20 jobs = 10 replicates for the Control network + 10 for the MS network.
# Each job is single-threaded (the network construction in the R package does
# not parallelise effectively), so they pack onto one node efficiently.
#
#   1-10   MSvC_Control
#   11-20  MSvC_MS
#
# Wall-clock estimate: all 20 jobs run concurrently, so the total is roughly
# the time of a single network build (~1.5-2 h at 800 genes x 10 networks x
# 500 cells on a modern core), plus the knockout step (~10 min).
# =============================================================================
#SBATCH --job-name=vko-stab
#SBATCH --array=1-20
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/stability_%A_%a.log
#SBATCH --error=logs/stability_%A_%a.err

set -euo pipefail
mkdir -p logs

if [ "${SLURM_ARRAY_TASK_ID}" -le 10 ]; then
  TAG="MSvC_Control"
  REP="${SLURM_ARRAY_TASK_ID}"
else
  TAG="MSvC_MS"
  REP=$(( SLURM_ARRAY_TASK_ID - 10 ))
fi

echo "tag=${TAG} rep=${REP} host=$(hostname) start=$(date -Is)"
Rscript code/03_controls/28_stability_replicates.R "${TAG}" "${REP}"
echo "done=$(date -Is)"

# Submit (the two conditions could also be separate submissions):
#   mkdir -p logs && sbatch code/03_controls/run_stability_slurm.sh
#
# Throughput tip: if the queue is long, submit with a core request of 1 and
# allow the scheduler to interleave the array; each job writes its own file and
# skips work that already exists, so the script is safe to re-run.
