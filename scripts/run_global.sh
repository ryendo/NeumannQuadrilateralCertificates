#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
: "${INTLAB_ROOT:?Set INTLAB_ROOT to the INTLAB installation directory.}"
: "${VEIGS_ROOT:?Set VEIGS_ROOT to the pinned veigs checkout.}"
mkdir -p "${ROOT}/results/global"
matlab -nodisplay -batch "addpath(genpath('${INTLAB_ROOT}')); cd '${ROOT}'; startintlab; r=QuadrilateralProofRunner('${INTLAB_ROOT}','${VEIGS_ROOT}'); r.setup(); r.runGlobal()" \
  | tee "${ROOT}/results/global/run.log"
