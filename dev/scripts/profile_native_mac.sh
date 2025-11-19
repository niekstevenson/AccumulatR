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

PROFILE_FILE=${UUBER_PROFILE_FILE:-/tmp/uuber_profile_r.txt}
PROFILE_DURATION=${UUBER_PROFILE_DURATION:-20}
PROFILE_INTERVAL=${UUBER_PROFILE_INTERVAL:-1}

export UUBER_PROFILE_TRIALS=${UUBER_PROFILE_TRIALS:-2000}
export UUBER_PROFILE_SEED=${UUBER_PROFILE_SEED:-2025}
export UUBER_REPO_ROOT="$REPO_ROOT"

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

sleep 1

if ! ps -p $R_PID >/dev/null 2>&1; then
  echo "R process exited before sampling could start" >&2
  wait $R_PID || true
  exit 1
fi

echo "Sampling PID $R_PID for ${PROFILE_DURATION}s (interval ${PROFILE_INTERVAL}ms) -> ${PROFILE_FILE}" >&2
if ! sample $R_PID $PROFILE_DURATION $PROFILE_INTERVAL -file "$PROFILE_FILE"; then
  echo "Warning: 'sample' exited with a non-zero status" >&2
fi

wait $R_PID || true

trap - EXIT

echo "Profiling complete. Output saved to ${PROFILE_FILE}"
