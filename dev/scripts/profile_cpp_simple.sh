#!/usr/bin/env bash

set -euo pipefail

if ! command -v sample >/dev/null 2>&1; then
  echo "Error: macOS 'sample' tool not found in PATH" >&2
  exit 1
fi

if ! command -v R >/dev/null 2>&1; then
  echo "Error: R not found in PATH" >&2
  exit 1
fi

SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd -P)
REPO_ROOT=$(cd -- "${SCRIPT_DIR}/../.." && pwd -P)

PROFILE_R_SCRIPT="${REPO_ROOT}/dev/scripts/profile.R"
if [[ ! -f "$PROFILE_R_SCRIPT" ]]; then
  echo "Error: profiling workload script not found at $PROFILE_R_SCRIPT" >&2
  exit 1
fi

PROFILE_FILE=${ACCUMULATR_PROFILE_FILE:-"$REPO_ROOT/dev/accumulatr_profile_cpp.txt"}
PROFILE_DURATION=${ACCUMULATR_PROFILE_DURATION:-20}
PROFILE_INTERVAL=${ACCUMULATR_PROFILE_INTERVAL:-1}

PROFILE_LIB=$(mktemp -d /tmp/accumulatr_profile_lib.XXXXXX)
PKG_CXXFLAGS="-O2 -g -fno-omit-frame-pointer" \
PKG_CFLAGS="-O2 -g -fno-omit-frame-pointer" \
PKG_CPPFLAGS="-g" \
R CMD INSTALL --preclean --no-multiarch --with-keep.source --library="${PROFILE_LIB}" "${REPO_ROOT}"
export R_LIBS_USER="${PROFILE_LIB}${R_LIBS_USER:+":${R_LIBS_USER}"}"

cd "$REPO_ROOT"

echo "Launching R profiling workload (profile.R)..."
R --vanilla --quiet -f "$PROFILE_R_SCRIPT" &
R_PID=$!

cleanup() {
  if ps -p $R_PID >/dev/null 2>&1; then
    kill $R_PID >/dev/null 2>&1 || true
  fi
  rm -rf "$PROFILE_LIB"
}
trap cleanup EXIT

if ! ps -p $R_PID >/dev/null 2>&1; then
  echo "R process exited before sampling could start" >&2
  wait $R_PID || true
  exit 1
fi

echo "Sampling PID $R_PID for ${PROFILE_DURATION}s (interval ${PROFILE_INTERVAL}ms) -> ${PROFILE_FILE}" >&2
sample_args=("$R_PID" "$PROFILE_DURATION")
if [[ "${PROFILE_INTERVAL}" -gt 0 ]]; then
  sample_args+=("$PROFILE_INTERVAL")
fi
sample_ok=0
for attempt in 1 2 3; do
  if sample "${sample_args[@]}" -file "$PROFILE_FILE"; then
    sample_ok=1
    break
  fi
  sleep 0.2
done
if [[ $sample_ok -ne 1 ]]; then
  echo "Warning: 'sample' exited with a non-zero status" >&2
fi

wait $R_PID || true

trap - EXIT

echo "Profiling complete. Output saved to ${PROFILE_FILE}"
