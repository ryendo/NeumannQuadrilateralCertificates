#!/usr/bin/env bash
set -euo pipefail

WORKERS="${1:-16}"
FACE_SUBDIVISIONS="${2:-12}"
OUTDIR="${3:-results/local_new}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"

if [[ -z "${INTLAB_ROOT:-}" && -z "${INTLAB_ROOT_PATTERN:-}" ]]; then
  echo "Set INTLAB_ROOT or INTLAB_ROOT_PATTERN." >&2
  exit 2
fi
mkdir -p "${ROOT}/${OUTDIR}"
pids=()
for k in $(seq 1 "${WORKERS}"); do
  if [[ -n "${INTLAB_ROOT_PATTERN:-}" ]]; then
    intlab_root="$(printf "${INTLAB_ROOT_PATTERN}" "$(( (k-1)%60+1 ))")"
  else
    intlab_root="${INTLAB_ROOT}"
  fi
  "${ROOT}/scripts/matlab_runner.sh" "${intlab_root}" \
    "maxNumCompThreads(1); qn_local_certificate_cover(${k},${WORKERS},${FACE_SUBDIVISIONS},fullfile('${ROOT}','${OUTDIR}'))" \
    > "${ROOT}/${OUTDIR}/log_${k}.txt" 2>&1 &
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
echo "Completed ${WORKERS} local-certificate workers."
