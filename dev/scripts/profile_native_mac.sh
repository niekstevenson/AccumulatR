#!/usr/bin/env bash

set -euo pipefail

if ! command -v sample >/dev/null 2>&1; then
  echo "Error: macOS 'sample' tool not found in PATH" >&2
  exit 1
fi

SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd -P)
REPO_ROOT=$(cd -- "${SCRIPT_DIR}/../.." && pwd -P)

PROFILE_R_SCRIPT="${REPO_ROOT}/dev/scripts/profile_examples.R"
if [[ ! -f "$PROFILE_R_SCRIPT" ]]; then
  echo "Error: profiling workload script not found at $PROFILE_R_SCRIPT" >&2
  exit 1
fi

PROFILE_FILE=${UUBER_PROFILE_FILE:-"$REPO_ROOT/dev/uuber_profile_r.txt"}
PROFILE_DURATION=${UUBER_PROFILE_DURATION:-20}
PROFILE_INTERVAL=${UUBER_PROFILE_INTERVAL:-1}

export UUBER_PROFILE_TRIALS=${UUBER_PROFILE_TRIALS:-2000}
export UUBER_PROFILE_SEED=${UUBER_PROFILE_SEED:-2025}
export UUBER_REPO_ROOT="$REPO_ROOT"
export UUBER_SKIP_DEVTOOLS=1

# Rebuild AccumulatR with profiling-friendly flags into a temp library
PROFILE_LIB=$(mktemp -d /tmp/uuber_profile_lib.XXXXXX)
PKG_CXXFLAGS="-O2 -g -fno-omit-frame-pointer" \
PKG_CFLAGS="-O2 -g -fno-omit-frame-pointer" \
PKG_CPPFLAGS="-g" \
R CMD INSTALL --preclean --no-multiarch --with-keep.source --library="${PROFILE_LIB}" "${REPO_ROOT}"
export R_LIBS_USER="${PROFILE_LIB}${R_LIBS_USER:+":${R_LIBS_USER}"}"

cd "$REPO_ROOT"

echo "Launching R profiling workload (examples 1-3)..."
R --vanilla --quiet -f "$PROFILE_R_SCRIPT" &
R_PID=$!

cleanup() {
  if ps -p $R_PID >/dev/null 2>&1; then
    kill $R_PID >/dev/null 2>&1 || true
  fi
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
  # Give the R process a moment to finish starting up before retrying
  sleep 0.2
done
if [[ $sample_ok -ne 1 ]]; then
  echo "Warning: 'sample' exited with a non-zero status" >&2
fi

wait $R_PID || true

trap - EXIT

echo "Profiling complete. Output saved to ${PROFILE_FILE}"