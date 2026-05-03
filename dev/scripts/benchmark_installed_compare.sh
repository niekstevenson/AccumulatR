#!/usr/bin/env bash

set -euo pipefail

if ! command -v R >/dev/null 2>&1; then
  echo "Error: R not found in PATH" >&2
  exit 1
fi

if ! command -v Rscript >/dev/null 2>&1; then
  echo "Error: Rscript not found in PATH" >&2
  exit 1
fi

if ! command -v git >/dev/null 2>&1; then
  echo "Error: git not found in PATH" >&2
  exit 1
fi

if ! command -v rsync >/dev/null 2>&1; then
  echo "Error: rsync not found in PATH" >&2
  exit 1
fi

SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd -P)
REPO_ROOT=$(cd -- "${SCRIPT_DIR}/../.." && pwd -P)
OUT_DIR="${ACCUMULATR_COMPARE_OUT_DIR:-"${REPO_ROOT}/dev/scripts/scratch_outputs"}"
BASELINE_REF="${ACCUMULATR_COMPARE_BASELINE_REF:-HEAD}"
TRIALS="${ACCUMULATR_COMPARE_TRIALS:-50 400}"
TARGET_SEC="${ACCUMULATR_COMPARE_TARGET_SEC:-0.12}"
MAX_INNER_REPS="${ACCUMULATR_COMPARE_MAX_INNER_REPS:-200}"
TIMESTAMP=$(date -u +"%Y%m%dT%H%M%SZ")
RUN_NAME="${ACCUMULATR_COMPARE_RUN_NAME:-installed_compare_${TIMESTAMP}}"

mkdir -p "$OUT_DIR"

WORK_ROOT=$(mktemp -d /tmp/acc_installed_compare.XXXXXX)
CANDIDATE_SRC="${WORK_ROOT}/candidate"
BASELINE_SRC="${WORK_ROOT}/baseline"
CANDIDATE_LIB="${WORK_ROOT}/candidate_lib"
BASELINE_LIB="${WORK_ROOT}/baseline_lib"
LOG_DIR="${OUT_DIR}/${RUN_NAME}_logs"
mkdir -p "$CANDIDATE_LIB" "$BASELINE_LIB" "$LOG_DIR"

cleanup() {
  if [[ "${ACCUMULATR_COMPARE_KEEP_WORK:-0}" == "1" ]]; then
    echo "Keeping temporary work directory: ${WORK_ROOT}" >&2
    return
  fi
  git -C "$REPO_ROOT" worktree remove --force "$BASELINE_SRC" >/dev/null 2>&1 || true
  rm -rf "$WORK_ROOT"
}
trap cleanup EXIT

candidate_commit=$(git -C "$REPO_ROOT" rev-parse HEAD)
candidate_dirty_count=$(git -C "$REPO_ROOT" status --porcelain | wc -l | tr -d ' ')
baseline_commit=$(git -C "$REPO_ROOT" rev-parse "$BASELINE_REF")

cat > "${LOG_DIR}/manifest.txt" <<EOF
run_name=${RUN_NAME}
timestamp_utc=${TIMESTAMP}
repo_root=${REPO_ROOT}
baseline_ref=${BASELINE_REF}
baseline_commit=${baseline_commit}
candidate_commit=${candidate_commit}
candidate_dirty_count=${candidate_dirty_count}
work_root=${WORK_ROOT}
candidate_src=${CANDIDATE_SRC}
baseline_src=${BASELINE_SRC}
candidate_lib=${CANDIDATE_LIB}
baseline_lib=${BASELINE_LIB}
trials=${TRIALS}
target_sec=${TARGET_SEC}
max_inner_reps=${MAX_INNER_REPS}
EOF

git -C "$REPO_ROOT" status --short > "${LOG_DIR}/candidate_git_status.txt"

echo "Copying candidate dirty tree to ${CANDIDATE_SRC}" >&2
mkdir -p "$CANDIDATE_SRC"
rsync -a \
  --exclude=".git" \
  --exclude="dev/scripts/scratch_outputs" \
  "${REPO_ROOT}/" \
  "${CANDIDATE_SRC}/"

echo "Creating baseline worktree ${BASELINE_SRC} at ${BASELINE_REF}" >&2
git -C "$REPO_ROOT" worktree add --detach "$BASELINE_SRC" "$BASELINE_REF" \
  > "${LOG_DIR}/baseline_worktree.log" 2>&1

echo "Installing candidate package" >&2
R CMD INSTALL --preclean -l "$CANDIDATE_LIB" "$CANDIDATE_SRC" \
  > "${LOG_DIR}/candidate_install.log" 2>&1

echo "Installing baseline package" >&2
R CMD INSTALL --preclean -l "$BASELINE_LIB" "$BASELINE_SRC" \
  > "${LOG_DIR}/baseline_install.log" 2>&1

run_benchmark() {
  local label=$1
  local src=$2
  local lib=$3
  local n_trials=$4
  local n_rep=$5
  local out_file=$6

  echo "Running ${label}, n_trials=${n_trials}, n_rep=${n_rep}" >&2
  (
    cd "$src"
    R_LIBS="$lib" \
    R_LIBS_USER="$lib" \
    ACCUMULATR_BENCH_INSTALLED=true \
    ACCUMULATR_BENCH_TRIALS="$n_trials" \
    ACCUMULATR_BENCH_N_REP="$n_rep" \
    ACCUMULATR_BENCH_TARGET_SEC="$TARGET_SEC" \
    ACCUMULATR_BENCH_MAX_INNER_REPS="$MAX_INNER_REPS" \
    ACCUMULATR_BENCH_OUT="$out_file" \
    Rscript dev/scripts/benchmark_speed.R
  ) > "${LOG_DIR}/${label}_n${n_trials}.log" 2>&1
}

for n_trials in $TRIALS; do
  case "$n_trials" in
    50) n_rep="${ACCUMULATR_COMPARE_N_REP_50:-5}" ;;
    400) n_rep="${ACCUMULATR_COMPARE_N_REP_400:-3}" ;;
    *) n_rep="${ACCUMULATR_COMPARE_N_REP:-3}" ;;
  esac

  candidate_csv="${OUT_DIR}/${RUN_NAME}_candidate_n${n_trials}.csv"
  baseline_csv="${OUT_DIR}/${RUN_NAME}_baseline_n${n_trials}.csv"
  compare_csv="${OUT_DIR}/${RUN_NAME}_compare_n${n_trials}.csv"

  run_benchmark "candidate" "$CANDIDATE_SRC" "$CANDIDATE_LIB" \
    "$n_trials" "$n_rep" "$candidate_csv"
  run_benchmark "baseline" "$BASELINE_SRC" "$BASELINE_LIB" \
    "$n_trials" "$n_rep" "$baseline_csv"

  Rscript - "$baseline_csv" "$candidate_csv" "$compare_csv" "$n_trials" <<'RS'
args <- commandArgs(trailingOnly = TRUE)
baseline_file <- args[[1]]
candidate_file <- args[[2]]
compare_file <- args[[3]]
n_trials <- as.integer(args[[4]])

baseline <- utils::read.csv(baseline_file, stringsAsFactors = FALSE)
candidate <- utils::read.csv(candidate_file, stringsAsFactors = FALSE)

required <- c("label", "median_per_eval_ms", "inner_reps", "git_commit", "git_dirty")
missing_baseline <- setdiff(required, names(baseline))
missing_candidate <- setdiff(required, names(candidate))
if (length(missing_baseline) != 0L || length(missing_candidate) != 0L) {
  stop(
    "Benchmark CSV missing required columns: baseline=",
    paste(missing_baseline, collapse = ","),
    " candidate=",
    paste(missing_candidate, collapse = ",")
  )
}

out <- merge(
  baseline[, c("label", "median_per_eval_ms", "inner_reps", "git_commit", "git_dirty")],
  candidate[, c("label", "median_per_eval_ms", "inner_reps", "git_commit", "git_dirty")],
  by = "label",
  suffixes = c("_baseline", "_candidate")
)
out$n_trials <- n_trials
out$ratio_candidate_over_baseline <- ifelse(
  out$median_per_eval_ms_baseline > 0,
  out$median_per_eval_ms_candidate / out$median_per_eval_ms_baseline,
  NA_real_
)
out$speedup_baseline_over_candidate <- ifelse(
  out$median_per_eval_ms_candidate > 0,
  out$median_per_eval_ms_baseline / out$median_per_eval_ms_candidate,
  NA_real_
)
out <- out[order(out$ratio_candidate_over_baseline, out$label), ]
utils::write.csv(out, compare_file, row.names = FALSE)

print(
  out[, c(
    "label",
    "median_per_eval_ms_baseline",
    "median_per_eval_ms_candidate",
    "ratio_candidate_over_baseline",
    "speedup_baseline_over_candidate"
  )],
  row.names = FALSE,
  digits = 4
)
RS
done

echo "Wrote installed comparison outputs under ${OUT_DIR}/${RUN_NAME}_*.csv" >&2
echo "Wrote logs under ${LOG_DIR}" >&2
