#!/usr/bin/env bash

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "$0")" && pwd -P)
repo_root=$(cd -- "$script_dir/../.." && pwd -P)
profile_lib=$(mktemp -d /tmp/accumulatr-profile-lib.XXXXXX)
start_file=$(mktemp /tmp/accumulatr-profile-start.XXXXXX)
rm "$start_file"

cleanup() {
  if [[ -n "${r_pid:-}" ]] && kill -0 "$r_pid" 2>/dev/null; then
    kill "$r_pid"
  fi
  rm -rf "$profile_lib"
  rm -f "$start_file"
}
trap cleanup EXIT

PKG_CXXFLAGS="-O2 -g -fno-omit-frame-pointer" \
PKG_CFLAGS="-O2 -g -fno-omit-frame-pointer" \
R CMD INSTALL --preclean --library="$profile_lib" "$repo_root"

cd "$repo_root"
mkdir -p dev/scripts/scratch_outputs
ACCUMULATR_PROFILE_START_FILE="$start_file" \
R_LIBS_USER="$profile_lib" \
Rscript dev/scripts/profile_workload_mixed.R &
r_pid=$!

while [[ ! -f "$start_file" ]]; do
  kill -0 "$r_pid"
  sleep 0.01
done

sample "$r_pid" 20 1 -file dev/scripts/scratch_outputs/profile_cpp.txt
wait "$r_pid"
