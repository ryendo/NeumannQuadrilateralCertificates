#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
: "${INTLAB_ROOT:?Set INTLAB_ROOT to the INTLAB installation directory.}"
matlab -nodisplay -batch "addpath(genpath('${INTLAB_ROOT}')); cd '${ROOT}'; startintlab; r=QuadrilateralProofRunner('${INTLAB_ROOT}'); r.setup(); r.smokeTest()"
