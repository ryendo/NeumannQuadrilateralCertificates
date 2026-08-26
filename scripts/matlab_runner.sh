#!/usr/bin/env bash
set -euo pipefail

INTLAB_ROOT="${1:?INTLAB root is required}"
MATLAB_CODE="${2:?MATLAB command is required}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
: "${VEIGS_ROOT:?Set VEIGS_ROOT to the pinned veigs checkout.}"

matlab -nodisplay -batch \
  "cd '${ROOT}'; addpath(fullfile('${ROOT}','src')); qn_setup('${INTLAB_ROOT}','${VEIGS_ROOT}'); ${MATLAB_CODE}"
