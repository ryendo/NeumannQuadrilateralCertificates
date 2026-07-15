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
for k in $(seq 1 "${WORKERS}"); do
  if [[ -n "${INTLAB_ROOT_PATTERN:-}" ]]; then
    intlab_root="$(printf "${INTLAB_ROOT_PATTERN}" "$(( (k-1)%60+1 ))")"
  else
    intlab_root="${INTLAB_ROOT}"
  fi
  matlab -nodisplay -batch "maxNumCompThreads(1); addpath(genpath('${intlab_root}')); cd '${ROOT}'; startintlab; r=QuadrilateralProofRunner('${intlab_root}'); r.setup(); r.runLocalWorker(${k},${WORKERS},${FACE_SUBDIVISIONS},fullfile('${ROOT}','${OUTDIR}'))" \
    > "${ROOT}/${OUTDIR}/log_${k}.txt" 2>&1 &
  sleep 1
done
echo "Launched ${WORKERS} local-certificate workers."
