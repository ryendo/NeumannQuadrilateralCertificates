#!/usr/bin/env bash
# Shared provenance and fresh-output checks for certificate launchers.

qn_prepare_run() {
  local root="$1"
  local outdir="$2"
  local certificate="$3"
  local workers="$4"
  local root_real output_real source_commit source_dirty
  local expected_veigs_commit actual_veigs_commit veigs_dirty source_code_sha256

  : "${VEIGS_ROOT:?Set VEIGS_ROOT to the pinned veigs checkout.}"
  for checked_path in "${root}" "${outdir}" "${VEIGS_ROOT}" \
      "${INTLAB_ROOT:-}" "${INTLAB_ROOT_PATTERN:-}"; do
    if [[ "${checked_path}" == *"'"* ]]; then
      echo "Paths containing an apostrophe are unsupported: ${checked_path}" >&2
      return 2
    fi
  done

  root_real="$(realpath "${root}")"
  case "${outdir}" in
    results/*) ;;
    *) echo "Output directory must be below results/: ${outdir}" >&2; return 2 ;;
  esac
  output_real="$(realpath -m "${root_real}/${outdir}")"
  case "${output_real}/" in
    "${root_real}/"*) ;;
    *) echo "Output directory must be inside the repository: ${outdir}" >&2; return 2 ;;
  esac
  if [[ -e "${output_real}" ]]; then
    echo "Refusing to reuse an existing result directory: ${output_real}" >&2
    return 2
  fi

  source_commit="$(git -C "${root_real}" rev-parse HEAD)"
  if [[ -n "${QN_SOURCE_COMMIT:-}" && "${QN_SOURCE_COMMIT}" != "${source_commit}" ]]; then
    echo "QN_SOURCE_COMMIT does not match the executing checkout." >&2
    return 2
  fi
  source_dirty=0
  while IFS= read -r status_line; do
    [[ -z "${status_line}" ]] && continue
    [[ "${status_line}" == '?? results/'* ]] && continue
    source_dirty=1
    break
  done < <(git -C "${root_real}" status --porcelain --untracked-files=all)
  if [[ ! "${source_commit}" =~ ^[0-9a-f]{40}$ ]]; then
    echo "QN_SOURCE_COMMIT must be a full lowercase Git commit." >&2
    return 2
  fi
  if [[ "${source_dirty}" != 0 ]]; then
    echo "Certificate runs require a clean tracked source tree." >&2
    return 2
  fi

  expected_veigs_commit='6556d39a0d9819bb172d232062b698aa76e420f6'
  actual_veigs_commit="$(git -C "${VEIGS_ROOT}" rev-parse HEAD)"
  if [[ -n "$(git -C "${VEIGS_ROOT}" status --porcelain --untracked-files=all)" ]]; then
    veigs_dirty=1
  else
    veigs_dirty=0
  fi
  if [[ "${actual_veigs_commit}" != "${expected_veigs_commit}" || "${veigs_dirty}" != 0 ]]; then
    echo "veigs must be the clean pinned checkout ${expected_veigs_commit}." >&2
    return 2
  fi

  export QN_SOURCE_COMMIT="${source_commit}"
  export QN_SOURCE_DIRTY=0
  export QN_VEIGS_COMMIT="${actual_veigs_commit}"
  export QN_RUN_ID="${QN_RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)-${certificate}-$$}"
  export QN_OUTPUT_PATH="${output_real}"
  source_code_sha256="$(git -C "${root_real}" ls-tree -r --full-tree HEAD -- src scripts tests data \
    | LC_ALL=C sort | sha256sum | awk '{print $1}')"
  mkdir "${QN_OUTPUT_PATH}"
  {
    echo "schema_version=1"
    echo "certificate=${certificate}"
    echo "run_id=${QN_RUN_ID}"
    echo "source_commit=${QN_SOURCE_COMMIT}"
    echo "source_dirty=${QN_SOURCE_DIRTY}"
    echo "source_code_sha256=${source_code_sha256}"
    echo "veigs_commit=${QN_VEIGS_COMMIT}"
    echo "veigs_root=${VEIGS_ROOT}"
    if [[ -n "${INTLAB_ROOT_PATTERN:-}" ]]; then
      echo "intlab_root_pattern=${INTLAB_ROOT_PATTERN}"
    else
      echo "intlab_root=${INTLAB_ROOT}"
    fi
    echo "workers=${workers}"
    echo "matlab=$(command -v matlab)"
    echo "host=$(hostname)"
    echo "pbs_job_id=${PBS_JOBID:-none}"
    echo "started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  } > "${QN_OUTPUT_PATH}/RUN_PROVENANCE.txt"
}

qn_finalize_run() {
  local output_path="$1"
  local completed_provenance
  completed_provenance="${output_path}/.RUN_PROVENANCE.complete"
  cp "${output_path}/RUN_PROVENANCE.txt" "${completed_provenance}"
  {
    echo "ended_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "status=complete"
  } >> "${completed_provenance}"
  (
    cd "${output_path}"
    # Worker logs are diagnostic and may be large; the certificate manifest
    # authenticates the result, metadata, completion markers, and summary files.
    find . -maxdepth 1 -type f ! -name 'log_*.txt' \
      ! -name SHA256SUMS ! -name SHA256SUMS.tmp \
      ! -name RUN_PROVENANCE.txt ! -name .RUN_PROVENANCE.complete -print0 \
      | sort -z | xargs -0 sha256sum > SHA256SUMS.tmp
    sha256sum -c SHA256SUMS.tmp
    provenance_sha256="$(sha256sum .RUN_PROVENANCE.complete | awk '{print $1}')"
    mv .RUN_PROVENANCE.complete RUN_PROVENANCE.txt
    printf '%s  %s\n' "${provenance_sha256}" RUN_PROVENANCE.txt >> SHA256SUMS.tmp
    mv SHA256SUMS.tmp SHA256SUMS
    sha256sum -c SHA256SUMS
  )
}
