#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "${ROOT}/scripts/run_provenance.sh"
OUTDIR="${1:-results/global_new}"
: "${INTLAB_ROOT:?Set INTLAB_ROOT to the INTLAB installation directory.}"
qn_prepare_run "${ROOT}" "${OUTDIR}" global 1
{
  echo 'n_init=3'
  echo 'max_depth=60'
} >> "${QN_OUTPUT_PATH}/RUN_PROVENANCE.txt"
"${ROOT}/scripts/matlab_runner.sh" "${INTLAB_ROOT}" \
  "qn_global_certified_cover(3,60,true,fullfile('${QN_OUTPUT_PATH}','worker_001.json')); r=qn_merge_global_results('${QN_OUTPUT_PATH}',1,fullfile('${QN_OUTPUT_PATH}','summary.json')); assert(r.complete && r.current_provenance_valid);" \
  | tee "${QN_OUTPUT_PATH}/run.log"
qn_finalize_run "${QN_OUTPUT_PATH}"
