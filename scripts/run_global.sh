#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
: "${INTLAB_ROOT:?Set INTLAB_ROOT to the INTLAB installation directory.}"
mkdir -p "${ROOT}/results/global"
"${ROOT}/scripts/matlab_runner.sh" "${INTLAB_ROOT}" \
  "qn_global_certified_cover(3,60,true,fullfile('${ROOT}','results','global','summary.json'))" \
  | tee "${ROOT}/results/global/run.log"
