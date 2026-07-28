#!/usr/bin/env bash
set -euo pipefail

GLOBAL_RESULTS="${1:?global results directory is required}"
WORKERS="${2:-48}"
FACE_SUBDIVISIONS="${3:-12}"
OUTDIR="${4:-results/local_rigorous_new}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"

if [[ -z "${INTLAB_ROOT_PATTERN:-}" ]]; then
  echo "Set INTLAB_ROOT_PATTERN before launching this watcher." >&2
  exit 2
fi
: "${VEIGS_ROOT:?Set VEIGS_ROOT to the pinned veigs checkout.}"

count_markers() {
  local directory="$1" pattern="$2"
  find "$directory" -maxdepth 1 -name "$pattern" -print 2>/dev/null | wc -l | tr -d ' '
}

GLOBAL_MARKER_PATTERN="${GLOBAL_MARKER_PATTERN:-worker_*.json}"
echo "Waiting for ${WORKERS_GLOBAL:-10} global worker JSON files in ${GLOBAL_RESULTS}."
while [[ "$(count_markers "$GLOBAL_RESULTS" "$GLOBAL_MARKER_PATTERN")" -lt "${WORKERS_GLOBAL:-10}" ]]; do
  date
  sleep 60
done

mkdir -p "$ROOT/$OUTDIR"
echo "Global workers finished; launching the hardened local certificate."
"$ROOT/scripts/run_local_workers.sh" "$WORKERS" "$FACE_SUBDIVISIONS" "$OUTDIR"

while [[ "$(count_markers "$ROOT/$OUTDIR" 'done_*.txt')" -lt "$WORKERS" ]]; do
  date
  sleep 60
done

SUMMARY_INTLAB_ROOT="$(printf "$INTLAB_ROOT_PATTERN" 1)"
matlab -nodisplay -batch "maxNumCompThreads(1); addpath(genpath('${SUMMARY_INTLAB_ROOT}')); startintlab; addpath('${ROOT}/src/local'); s=qn_summarize_local_results('${ROOT}/${OUTDIR}'); save('${ROOT}/${OUTDIR}/summary.mat','s'); disp(s);" \
  > "$ROOT/$OUTDIR/summary.log" 2>&1
echo "Local certificate and summary completed."
