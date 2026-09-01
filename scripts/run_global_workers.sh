#!/usr/bin/env bash
set -euo pipefail

WORKERS="${1:-10}"
OUTDIR="${2:-results/global_new}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
: "${INTLAB_ROOT_PATTERN:?Set INTLAB_ROOT_PATTERN, e.g. /path/Intlab_V12_no%d}"
mkdir -p "${ROOT}/${OUTDIR}"

pids=()
for k in $(seq 1 "${WORKERS}"); do
  intlab_root="$(printf "${INTLAB_ROOT_PATTERN}" "$k")"
  "${ROOT}/scripts/matlab_runner.sh" "${intlab_root}" \
    "maxNumCompThreads(1); qn_global_certified_cover(3,60,true,fullfile('${ROOT}','${OUTDIR}',sprintf('worker_%03d.json',${k})),${k},${WORKERS})" \
    > "${ROOT}/${OUTDIR}/log_${k}.txt" 2>&1 &
  pids+=("$!")
  sleep 1
done
failed=0
for pid in "${pids[@]}"; do
  if ! wait "${pid}"; then failed=1; fi
done
if (( failed != 0 )); then
  echo "At least one global-certificate worker failed; inspect ${ROOT}/${OUTDIR}/log_*.txt." >&2
  exit 1
fi
echo "Completed ${WORKERS} workers for the initial boxes in Appendix B."
echo "After completion, merge with qn_merge_global_results('${OUTDIR}',${WORKERS},'${OUTDIR}/summary.json')."
