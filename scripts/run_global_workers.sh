#!/usr/bin/env bash
set -euo pipefail

WORKERS="${1:-10}"
OUTDIR="${2:-results/global_new}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "${ROOT}/scripts/run_provenance.sh"
: "${INTLAB_ROOT_PATTERN:?Set INTLAB_ROOT_PATTERN, e.g. /path/Intlab_V12_no%d}"
qn_prepare_run "${ROOT}" "${OUTDIR}" global "${WORKERS}"
{
  echo 'n_init=3'
  echo 'max_depth=60'
} >> "${QN_OUTPUT_PATH}/RUN_PROVENANCE.txt"

pids=()
for k in $(seq 1 "${WORKERS}"); do
  intlab_root="$(printf "${INTLAB_ROOT_PATTERN}" "$k")"
  "${ROOT}/scripts/matlab_runner.sh" "${intlab_root}" \
    "maxNumCompThreads(1); qn_global_certified_cover(3,60,true,fullfile('${QN_OUTPUT_PATH}',sprintf('worker_%03d.json',${k})),${k},${WORKERS})" \
    > "${QN_OUTPUT_PATH}/log_${k}.txt" 2>&1 &
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
summary_intlab_root="$(printf "${INTLAB_ROOT_PATTERN}" 1)"
"${ROOT}/scripts/matlab_runner.sh" "${summary_intlab_root}" \
  "r=qn_merge_global_results('${QN_OUTPUT_PATH}',${WORKERS},fullfile('${QN_OUTPUT_PATH}','summary.json')); assert(r.complete && r.current_provenance_valid);"
qn_finalize_run "${QN_OUTPUT_PATH}"
echo "Completed ${WORKERS} workers for the initial boxes in Appendix B."
