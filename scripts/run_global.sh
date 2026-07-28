#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
: "${INTLAB_ROOT:?Set INTLAB_ROOT to the INTLAB installation directory.}"
mkdir -p "${ROOT}/results/global"
"${ROOT}/scripts/matlab_runner.sh" "${INTLAB_ROOT}" "r.runGlobal()" \
  | tee "${ROOT}/results/global/run.log"
