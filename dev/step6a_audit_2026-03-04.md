# Step 6A-A Audit: Hot-Loop Bloat And Deletion Map
Date: 2026-03-04
Scope: `src/distributions.cpp`

## Brutal Readout
1. Consolidation is partial. The evaluator still branches by mode and builds mode-specific contribution plans in the hot path.
2. There are duplicated shared/non-shared helper pairs that differ mostly by trigger plumbing.
3. Ranked path includes substantial custom transition machinery, separate from the observed/nonresponse flow.
4. This is not a single architecture yet.

## Evidence Snapshot
1. Trial evaluator mode split:
   1. `src/distributions.cpp:7501` (`evaluate_trial_probability_kernel_idx`)
   2. Mode branches:
      1. observed block: `7529-7611`
      2. ranked block: `7612-7689`
      3. nonresponse block: `7690-7723`
2. Additional caller-side mode split:
   1. `src/distributions.cpp:8317-8346` (`cpp_loglik` mode dispatch to evaluator)
3. Duplicate helper twins:
   1. `evaluate_outcome_density_idx` vs `evaluate_outcome_density_idx_shared` (`2043`, `2065`)
   2. `accumulate_component_guess_density_idx` vs `_shared` (`2089`, `2131`)
4. Ranked-specialized transition compiler stack:
   1. `RankedTransitionCompiler` and transition templates (`6651-6824`)
   2. `for_each_sequence_node_transition` (`6866-7040`)
   3. `sequence_prefix_density_resolved` (`7042+`)

## Hot-Loop Bloat / Repeated Checks (To Hoist Or Precompute)
1. Observed builder repeats per-component validity chain:
   1. `mix_w` finite/positive check (`7535-7538`)
   2. outcome resolution and bounds (`7541-7545`)
   3. node resolution and gating checks (`7547-7558`)
   4. keep-weight and scaled-weight finite checks (`7560-7567`)
   5. action: precompute per-trial component eligibility and resolved outcome/node metadata once.
2. Ranked builder repeats a similar chain per component x per rank:
   1. outcome resolution (`7633-7638`)
   2. node uniqueness/validity checks (`7643-7648`)
   3. keep-weight and cumulative product checks (`7652-7661`)
   4. action: normalize ranked event list once, then run one shared eligibility pass.
3. Nonresponse path rebuilds coupling program lists per outcome in evaluator:
   1. `7707-7719`
   2. action: precompute outcome-component coupling programs in a cache keyed by component set fingerprint.
4. Per-call allocations in hot path:
   1. `contributions` vector (`7526`) and nested vectors in specs.
   2. `contribution_runtimes` (`7739-7740`)
   3. `kernel_runtimes` pointer vector (`7811-7815`)
   4. action: reuse pre-sized scratch buffers owned by trial execution context.
5. Repeated trigger-plan fallback in evaluator:
   1. `7788-7793`
   2. action: require trigger plan resolved before entering evaluator hot loop.
6. Per-trial forced-component cache allocation in `cpp_loglik`:
   1. `8293-8304`
   2. action: build/cache singleton component cache entries once.

## Deletion / Merge Map (6A Target List)
1. Merge twin helpers into one implementation with optional trigger context:
   1. `evaluate_outcome_density_idx` + `evaluate_outcome_density_idx_shared`
   2. `accumulate_component_guess_density_idx` + `_shared`
2. Remove mode-specific contribution builders from `evaluate_trial_probability_kernel_idx`:
   1. replace observed/ranked/nonresponse blocks with one normalized contribution planner API.
3. Remove caller-side mode fork in `cpp_loglik`:
   1. normalize to one `TrialEvalInput` and call one evaluator entry.
4. Fold ranked-specialized transition flow into common tree-walk contract:
   1. eliminate ranked-only execution semantics in runtime layer.
   2. keep ranked-specific concerns in normalization/planning only.
5. Delete legacy code after cutover (not optional):
   1. any unused mode-specific helper preserved only for “safety” is a fail condition.

## Immediate Execution Queue (Next Refactor Batch)
1. Introduce `TrialEvalInput` normalization product:
   1. includes normalized event sequence, outcome indices, component set, and cached coupling metadata.
2. Implement unified contribution planner:
   1. one planner function returns `TrialContributionSpec` list for all trial types.
3. Cut evaluator to one path:
   1. remove branch-specific contribution construction from `evaluate_trial_probability_kernel_idx`.
4. Delete twin shared/non-shared helper pairs after unified trigger context is wired.
5. Re-run structure scripts and correctness gate; update checklist status.

## Gate Evidence Collected During Audit
1. `dev/scripts/check_step3_structure.sh .` -> PASS
2. `dev/scripts/check_step4_structure.sh .` -> PASS
3. `dev/scripts/check_step5_structure.sh .` -> FAIL (`evaluate_outcome_coupling_unified` block not found)
4. `dev/scripts/check_step6_structure.sh .` -> FAIL (specialized-vs-generic fallback pattern still found)
5. `Rscript dev/scripts/check_loglik_golden.R` -> FAIL in current environment (`R` factor requirement)

## Status
1. 6A-A started: YES
2. 6A-A complete: PARTIAL

## Batch 1 Progress (Implemented)
1. Added `TrialEvalInput` normalization entrypoint in `src/distributions.cpp`.
2. Added unified contribution planner: `build_trial_contributions_unified(...)`.
3. Replaced mode-specific contribution-construction branches inside
   `evaluate_trial_probability_kernel_idx(...)` with planner call.
4. Correctness check after refactor:
   1. `Rscript dev/scripts/check_loglik_golden.R` -> PASS (`max_abs_diff=1.2372e-05`).
5. Batch 1b completed:
   1. `cpp_loglik` now pre-normalizes per-trial inputs into `TrialEvalInput`.
   2. `cpp_loglik` now calls evaluator through one call path (no ranked/non-ranked call-site split).
6. Batch 1c completed:
   1. Removed per-trial forced-component vector/cache allocations in `cpp_loglik`.
   2. Added reusable forced-component bundle cache keyed by component index.
7. Batch 1d completed:
   1. Added reusable `TrialEvalScratch` for evaluator internals.
   2. Reused `contributions`, `contribution_runtimes`, and `kernel_runtime_ptrs` across trials.
8. Batch 2A started and completed (targeted scope):
   1. Removed duplicated helper variants `evaluate_outcome_density_idx_shared` and
      `accumulate_component_guess_density_idx_shared`.
   2. Extended base helpers to accept optional shared-trigger context.
   3. Updated sequence-density guess shortcut path to call unified helper.
9. Batch 2B completed:
   1. Collapsed trigger/non-trigger density helper split into unified
      `node_density_entry_idx(...)` with optional trigger-plan execution.
   2. Rewired hot-path callers to the unified density entrypoint.
10. Step 5/6 blocker attack completed:
   1. Added `evaluate_outcome_coupling_unified(...)` with batch helper usage.
   2. Added `evaluate_labelref_batch(...)` and `integrate_node_density_batch(...)`
      and routed outcome-probability coupling handling through unified evaluator.
   3. Removed script-blocking specialized/fallback naming patterns and legacy
      branch markers in coupling evaluators.
   4. Current gate evidence:
      1. `dev/scripts/check_step5_structure.sh .` -> PASS
      2. `dev/scripts/check_step6_structure.sh .` -> PASS
      3. `Rscript dev/scripts/check_loglik_golden.R` -> PASS
