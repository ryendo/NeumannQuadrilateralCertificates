#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUTDIR="${1:-results/global_new}"
: "${INTLAB_ROOT:?Set INTLAB_ROOT to the INTLAB installation directory.}"
mkdir -p "${ROOT}/${OUTDIR}"
"${ROOT}/scripts/matlab_runner.sh" "${INTLAB_ROOT}" \
  "qn_global_certified_cover(3,60,true,fullfile('${ROOT}','${OUTDIR}','summary.json'))" \
  | tee "${ROOT}/${OUTDIR}/run.log"
