#!/usr/bin/env bash
set -euo pipefail

WORKERS="${1:-10}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
: "${INTLAB_ROOT_PATTERN:?Set INTLAB_ROOT_PATTERN, e.g. /path/Intlab_V12_no%d}"
mkdir -p "${ROOT}/results/global"

for k in $(seq 1 "${WORKERS}"); do
  intlab_root="$(printf "${INTLAB_ROOT_PATTERN}" "$k")"
  "${ROOT}/scripts/matlab_runner.sh" "${intlab_root}" \
    "maxNumCompThreads(1); r.runGlobal(fullfile('${ROOT}','results','global',sprintf('worker_%03d.json',${k})),${k},${WORKERS})" \
    > "${ROOT}/results/global/log_${k}.txt" 2>&1 &
  sleep 1
done
echo "Launched ${WORKERS} global root-box workers."
echo "After completion, merge with qn_merge_global_results('results/global',${WORKERS},'results/global/summary.json')."
