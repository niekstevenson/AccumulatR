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

PROFILE_R_SCRIPT=${UUBER_PROFILE_R_SCRIPT:-"${REPO_ROOT}/dev/scripts/profile_workload_nested.R"}
if [[ ! -f "$PROFILE_R_SCRIPT" ]]; then
  echo "Error: profiling workload script not found at $PROFILE_R_SCRIPT" >&2
  exit 1
fi

PROFILE_FILE=${UUBER_PROFILE_FILE:-"$REPO_ROOT/dev/uuber_profile_cpp.txt"}
PROFILE_DURATION=${UUBER_PROFILE_DURATION:-30}
PROFILE_INTERVAL=${UUBER_PROFILE_INTERVAL:-1}
PROFILE_START_TIMEOUT=${UUBER_PROFILE_START_TIMEOUT:-60}

PROFILE_TMP=$(mktemp -d /tmp/uuber_profile_run.XXXXXX)
PROFILE_START_FILE="${PROFILE_TMP}/profile.start"
PROFILE_END_FILE="${PROFILE_TMP}/profile.end"
PROFILE_WRAPPER_SCRIPT="${PROFILE_TMP}/profile_wrapper.R"

cat > "$PROFILE_WRAPPER_SCRIPT" <<'EOF'
start_file <- Sys.getenv("UUBER_PROFILE_START_FILE", "")
end_file <- Sys.getenv("UUBER_PROFILE_END_FILE", "")
target_script <- Sys.getenv("UUBER_PROFILE_R_SCRIPT", "")

if (!nzchar(target_script)) {
  stop("UUBER_PROFILE_R_SCRIPT not set")
}
if (nzchar(start_file)) {
  writeLines("start", start_file)
}
source(target_script)
if (nzchar(end_file)) {
  writeLines("end", end_file)
}
EOF

SAMPLE_PID=""
echo "Launching R profiling workload (${PROFILE_R_SCRIPT})..."
UUBER_PROFILE_START_FILE="$PROFILE_START_FILE" \
UUBER_PROFILE_END_FILE="$PROFILE_END_FILE" \
UUBER_PROFILE_R_SCRIPT="$PROFILE_R_SCRIPT" \
R --vanilla --quiet -f "$PROFILE_WRAPPER_SCRIPT" &
R_PID=$!

cleanup() {
  if [[ -n "$SAMPLE_PID" ]] && ps -p $SAMPLE_PID >/dev/null 2>&1; then
    kill -INT $SAMPLE_PID >/dev/null 2>&1 || true
    wait $SAMPLE_PID >/dev/null 2>&1 || true
  fi
  if ps -p $R_PID >/dev/null 2>&1; then
    kill $R_PID >/dev/null 2>&1 || true
  fi
  rm -rf "$PROFILE_TMP"
}
trap cleanup EXIT

if ! ps -p $R_PID >/dev/null 2>&1; then
  echo "R process exited before sampling could start" >&2
  wait $R_PID || true
  exit 1
fi

start_deadline=$((SECONDS + PROFILE_START_TIMEOUT))
echo "Waiting for workload start marker (timeout ${PROFILE_START_TIMEOUT}s)..." >&2
while [[ ! -f "$PROFILE_START_FILE" ]]; do
  if ! ps -p $R_PID >/dev/null 2>&1; then
    echo "R process exited before workload start marker" >&2
    wait $R_PID || true
    exit 1
  fi
  if [[ $SECONDS -ge $start_deadline ]]; then
    echo "Timed out waiting for workload start marker" >&2
    kill $R_PID >/dev/null 2>&1 || true
    wait $R_PID || true
    exit 1
  fi
  sleep 0.05
done

echo "Sampling PID $R_PID for ${PROFILE_DURATION}s (interval ${PROFILE_INTERVAL}ms) -> ${PROFILE_FILE}" >&2
sample_args=("$R_PID" "$PROFILE_DURATION")
if [[ "${PROFILE_INTERVAL}" -gt 0 ]]; then
  sample_args+=("$PROFILE_INTERVAL")
fi
sample "${sample_args[@]}" -file "$PROFILE_FILE" &
SAMPLE_PID=$!

while [[ ! -f "$PROFILE_END_FILE" ]]; do
  if ! ps -p $R_PID >/dev/null 2>&1; then
    break
  fi
  if ! ps -p $SAMPLE_PID >/dev/null 2>&1; then
    break
  fi
  sleep 0.05
done

if ps -p $SAMPLE_PID >/dev/null 2>&1; then
  kill -INT $SAMPLE_PID >/dev/null 2>&1 || true
fi
wait $SAMPLE_PID >/dev/null 2>&1 || true
wait $R_PID || true

trap - EXIT
rm -rf "$PROFILE_TMP"

echo "Profiling complete. Output saved to ${PROFILE_FILE}"
