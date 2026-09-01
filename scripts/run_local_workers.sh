#!/usr/bin/env bash
set -euo pipefail

WORKERS="${1:-16}"
FACE_SUBDIVISIONS="${2:-12}"
OUTDIR="${3:-results/local_new}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
source "${ROOT}/scripts/run_provenance.sh"

if [[ -z "${INTLAB_ROOT:-}" && -z "${INTLAB_ROOT_PATTERN:-}" ]]; then
  echo "Set INTLAB_ROOT or INTLAB_ROOT_PATTERN." >&2
  exit 2
fi
qn_prepare_run "${ROOT}" "${OUTDIR}" local "${WORKERS}"
{
  echo "face_subdivisions=${FACE_SUBDIVISIONS}"
  echo "expected_top_level_boxes=$((8*FACE_SUBDIVISIONS*FACE_SUBDIVISIONS*FACE_SUBDIVISIONS))"
} >> "${QN_OUTPUT_PATH}/RUN_PROVENANCE.txt"
pids=()
for k in $(seq 1 "${WORKERS}"); do
  if [[ -n "${INTLAB_ROOT_PATTERN:-}" ]]; then
    intlab_root="$(printf "${INTLAB_ROOT_PATTERN}" "$(( (k-1)%60+1 ))")"
  else
    intlab_root="${INTLAB_ROOT}"
  fi
  "${ROOT}/scripts/matlab_runner.sh" "${intlab_root}" \
    "maxNumCompThreads(1); qn_local_certificate_cover(${k},${WORKERS},${FACE_SUBDIVISIONS},'${QN_OUTPUT_PATH}')" \
    > "${QN_OUTPUT_PATH}/log_${k}.txt" 2>&1 &
  pids+=("$!")
  sleep 1
done
failed=0
for pid in "${pids[@]}"; do
  if ! wait "${pid}"; then failed=1; fi
done
if (( failed != 0 )); then
  echo "At least one local-certificate worker failed; inspect ${ROOT}/${OUTDIR}/log_*.txt." >&2
  exit 1
fi
if [[ -n "${INTLAB_ROOT_PATTERN:-}" ]]; then
  summary_intlab_root="$(printf "${INTLAB_ROOT_PATTERN}" 1)"
else
  summary_intlab_root="${INTLAB_ROOT}"
fi
"${ROOT}/scripts/matlab_runner.sh" "${summary_intlab_root}" \
  "s=qn_summarize_local_results('${QN_OUTPUT_PATH}',${WORKERS},${FACE_SUBDIVISIONS}); assert(s.verified); f=fopen(fullfile('${QN_OUTPUT_PATH}','summary.json'),'w'); assert(f>=0); fprintf(f,'%s\n',jsonencode(s,'PrettyPrint',true)); fclose(f);"
qn_finalize_run "${QN_OUTPUT_PATH}"
echo "Completed ${WORKERS} local-certificate workers."
