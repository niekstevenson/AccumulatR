#include "kernel_executor.h"

#include <cmath>
#include <limits>
#include <vector>

#include "kernel_program.h"

namespace uuber {
namespace {

inline double clamp_probability(double x) {
  if (!std::isfinite(x)) {
    return 0.0;
  }
  if (x < 0.0) {
    return 0.0;
  }
  if (x > 1.0) {
    return 1.0;
  }
  return x;
}

inline double safe_density(double x) {
  if (!std::isfinite(x) || x < 0.0) {
    return 0.0;
  }
  return x;
}

inline KernelNodeValues default_guard_eval(const KernelNodeValues &ref,
                                           const KernelNodeValues &block,
                                           EvalNeed need) {
  KernelNodeValues out;
  out.density = needs_density(need) ? safe_density(ref.density * block.survival)
                                    : 0.0;
  if (needs_cdf(need)) {
    out.cdf = clamp_probability(ref.cdf * block.survival);
  }
  if (needs_survival(need)) {
    out.survival =
        needs_cdf(need) ? clamp_probability(1.0 - out.cdf) : block.survival;
  }
  return out;
}

inline EvalNeed mask_to_need(std::uint8_t mask) {
  return static_cast<EvalNeed>(mask);
}

void prepare_runtime_for_plan(const KernelProgram &program,
                              KernelRuntimeState &runtime,
                              const KernelQueryPlan &plan) {
  if (!runtime.initialized || runtime.program != &program ||
      runtime.slots.size() != program.ops.size()) {
    reset_kernel_runtime(program, runtime);
  }
  if (runtime.active_plan != &plan) {
    invalidate_kernel_runtime_from_slot(runtime, 0);
    runtime.active_plan = &plan;
  }
}

bool eval_kernel_slots_incremental(const KernelProgram &program,
                                   KernelRuntimeState &runtime,
                                   const KernelQueryPlan &plan,
                                   const KernelEventEvaluator &event_eval,
                                   const KernelGuardEvaluator &guard_eval) {
  if (!plan.valid || plan.max_slot < 0 ||
      plan.slot_need_masks.size() != program.ops.size()) {
    return false;
  }
  const int eval_slot_count = plan.max_slot + 1;
  for (int op_idx = std::max(0, runtime.computed_upto + 1);
       op_idx <= plan.max_slot;
       ++op_idx) {
    const EvalNeed slot_need = mask_to_need(
        plan.slot_need_masks[static_cast<std::size_t>(op_idx)]);
    if (!has_any_need(slot_need)) {
      runtime.slots[static_cast<std::size_t>(op_idx)] = KernelNodeValues{};
      continue;
    }
    const KernelOp &op = program.ops[static_cast<std::size_t>(op_idx)];
    KernelNodeValues out;
    switch (op.code) {
    case KernelOpCode::Event: {
      out = event_eval(op.event_idx, slot_need);
      out.density = safe_density(out.density);
      out.survival = clamp_probability(out.survival);
      out.cdf = clamp_probability(out.cdf);
      break;
    }
    case KernelOpCode::And: {
      if (op.child_begin < 0 || op.child_count <= 0 ||
          op.child_begin + op.child_count >
              static_cast<int>(program.children.size())) {
        out = KernelNodeValues{};
        break;
      }
      const int n = op.child_count;
      for (int i = 0; i < n; ++i) {
        int slot = program.children[static_cast<std::size_t>(op.child_begin + i)];
        if (slot < 0 || slot >= eval_slot_count) {
          runtime.child_primary[static_cast<std::size_t>(i)] = 0.0;
          runtime.child_density[static_cast<std::size_t>(i)] = 0.0;
          continue;
        }
        runtime.child_primary[static_cast<std::size_t>(i)] =
            clamp_probability(runtime.slots[static_cast<std::size_t>(slot)].cdf);
        runtime.child_density[static_cast<std::size_t>(i)] =
            safe_density(runtime.slots[static_cast<std::size_t>(slot)].density);
      }
      runtime.prefix[0] = 1.0;
      runtime.suffix[static_cast<std::size_t>(n)] = 1.0;
      for (int i = 0; i < n; ++i) {
        runtime.prefix[static_cast<std::size_t>(i + 1)] =
            runtime.prefix[static_cast<std::size_t>(i)] *
            runtime.child_primary[static_cast<std::size_t>(i)];
      }
      for (int i = n - 1; i >= 0; --i) {
        runtime.suffix[static_cast<std::size_t>(i)] =
            runtime.suffix[static_cast<std::size_t>(i + 1)] *
            runtime.child_primary[static_cast<std::size_t>(i)];
      }
      const double cdf = clamp_probability(runtime.prefix[static_cast<std::size_t>(n)]);
      if (needs_cdf(slot_need)) {
        out.cdf = cdf;
      }
      if (needs_survival(slot_need)) {
        out.survival = clamp_probability(1.0 - cdf);
      }
      if (needs_density(slot_need)) {
        double d = 0.0;
        for (int i = 0; i < n; ++i) {
          const double others =
              runtime.prefix[static_cast<std::size_t>(i)] *
              runtime.suffix[static_cast<std::size_t>(i + 1)];
          d += runtime.child_density[static_cast<std::size_t>(i)] *
               clamp_probability(others);
        }
        out.density = safe_density(d);
      }
      break;
    }
    case KernelOpCode::Or: {
      if (op.child_begin < 0 || op.child_count <= 0 ||
          op.child_begin + op.child_count >
              static_cast<int>(program.children.size())) {
        out = KernelNodeValues{};
        break;
      }
      const int n = op.child_count;
      for (int i = 0; i < n; ++i) {
        int slot = program.children[static_cast<std::size_t>(op.child_begin + i)];
        if (slot < 0 || slot >= eval_slot_count) {
          runtime.child_primary[static_cast<std::size_t>(i)] = 1.0;
          runtime.child_density[static_cast<std::size_t>(i)] = 0.0;
          continue;
        }
        runtime.child_primary[static_cast<std::size_t>(i)] =
            clamp_probability(
                runtime.slots[static_cast<std::size_t>(slot)].survival);
        runtime.child_density[static_cast<std::size_t>(i)] =
            safe_density(runtime.slots[static_cast<std::size_t>(slot)].density);
      }
      runtime.prefix[0] = 1.0;
      runtime.suffix[static_cast<std::size_t>(n)] = 1.0;
      for (int i = 0; i < n; ++i) {
        runtime.prefix[static_cast<std::size_t>(i + 1)] =
            runtime.prefix[static_cast<std::size_t>(i)] *
            runtime.child_primary[static_cast<std::size_t>(i)];
      }
      for (int i = n - 1; i >= 0; --i) {
        runtime.suffix[static_cast<std::size_t>(i)] =
            runtime.suffix[static_cast<std::size_t>(i + 1)] *
            runtime.child_primary[static_cast<std::size_t>(i)];
      }
      const double survival =
          clamp_probability(runtime.prefix[static_cast<std::size_t>(n)]);
      if (needs_survival(slot_need)) {
        out.survival = survival;
      }
      if (needs_cdf(slot_need)) {
        out.cdf = clamp_probability(1.0 - survival);
      }
      if (needs_density(slot_need)) {
        double d = 0.0;
        for (int i = 0; i < n; ++i) {
          const double others =
              runtime.prefix[static_cast<std::size_t>(i)] *
              runtime.suffix[static_cast<std::size_t>(i + 1)];
          d += runtime.child_density[static_cast<std::size_t>(i)] *
               clamp_probability(others);
        }
        out.density = safe_density(d);
      }
      break;
    }
    case KernelOpCode::Not: {
      if (op.child_begin < 0 || op.child_count <= 0 ||
          op.child_begin >= static_cast<int>(program.children.size())) {
        out = KernelNodeValues{};
        break;
      }
      int child_slot = program.children[static_cast<std::size_t>(op.child_begin)];
      if (child_slot < 0 || child_slot >= eval_slot_count) {
        out = KernelNodeValues{};
        break;
      }
      const KernelNodeValues &child =
          runtime.slots[static_cast<std::size_t>(child_slot)];
      const double child_cdf = clamp_probability(child.cdf);
      if (needs_cdf(slot_need)) {
        out.cdf = clamp_probability(1.0 - child_cdf);
      }
      if (needs_survival(slot_need)) {
        out.survival = child_cdf;
      }
      out.density = 0.0;
      break;
    }
    case KernelOpCode::Guard: {
      const KernelNodeValues zero{};
      const KernelNodeValues &ref =
          (op.ref_slot >= 0 && op.ref_slot < eval_slot_count)
              ? runtime.slots[static_cast<std::size_t>(op.ref_slot)]
              : zero;
      const KernelNodeValues &block =
          (op.blocker_slot >= 0 && op.blocker_slot < eval_slot_count)
              ? runtime.slots[static_cast<std::size_t>(op.blocker_slot)]
              : zero;
      if (guard_eval) {
        out = guard_eval(op, ref, block, slot_need);
      } else {
        out = default_guard_eval(ref, block, slot_need);
      }
      out.density = safe_density(out.density);
      out.survival = clamp_probability(out.survival);
      out.cdf = clamp_probability(out.cdf);
      break;
    }
    }
    if (op.out_slot >= 0 && op.out_slot < eval_slot_count) {
      runtime.slots[static_cast<std::size_t>(op.out_slot)] = out;
    }
  }

  runtime.computed_upto = std::max(runtime.computed_upto, plan.max_slot);
  return true;
}

} // namespace

void reset_kernel_runtime(const KernelProgram &program,
                          KernelRuntimeState &runtime) {
  runtime.program = &program;
  runtime.active_plan = nullptr;
  runtime.slots.resize(program.ops.size());
  const int max_children = std::max(0, program.max_child_count);
  runtime.child_primary.resize(static_cast<std::size_t>(max_children));
  runtime.child_density.resize(static_cast<std::size_t>(max_children));
  runtime.prefix.resize(static_cast<std::size_t>(max_children + 1));
  runtime.suffix.resize(static_cast<std::size_t>(max_children + 1));
  runtime.computed_upto = -1;
  runtime.initialized = true;
}

void invalidate_kernel_runtime_from_slot(KernelRuntimeState &runtime,
                                         int slot_begin) {
  if (!runtime.initialized) {
    return;
  }
  if (slot_begin <= 0) {
    runtime.computed_upto = -1;
    return;
  }
  if (slot_begin - 1 < runtime.computed_upto) {
    runtime.computed_upto = slot_begin - 1;
  }
}

bool eval_kernel_query_plan_incremental(
    const KernelProgram &program, KernelRuntimeState &runtime,
    const KernelQueryPlan &plan, const KernelEventEvaluator &event_eval,
    const KernelGuardEvaluator &guard_eval,
    std::vector<KernelNodeValues> &out_values) {
  out_values.clear();
  if (!program.valid || !plan.valid || !event_eval ||
      plan.target_slots.empty()) {
    return false;
  }
  prepare_runtime_for_plan(program, runtime, plan);
  if (plan.max_slot > runtime.computed_upto) {
    if (!eval_kernel_slots_incremental(program, runtime, plan, event_eval,
                                       guard_eval)) {
      return false;
    }
  }
  out_values.reserve(plan.target_slots.size());
  for (int slot : plan.target_slots) {
    if (slot < 0 || slot >= static_cast<int>(runtime.slots.size())) {
      return false;
    }
    out_values.push_back(runtime.slots[static_cast<std::size_t>(slot)]);
  }
  return true;
}

bool eval_kernel_single_query_plan_incremental(
    const KernelProgram &program, KernelRuntimeState &runtime,
    const KernelQueryPlan &plan, const KernelEventEvaluator &event_eval,
    const KernelGuardEvaluator &guard_eval, KernelNodeValues &out_values) {
  if (!program.valid || !plan.valid || !event_eval ||
      plan.target_slots.size() != 1u) {
    return false;
  }
  prepare_runtime_for_plan(program, runtime, plan);
  if (plan.max_slot > runtime.computed_upto) {
    if (!eval_kernel_slots_incremental(program, runtime, plan, event_eval,
                                       guard_eval)) {
      return false;
    }
  }
  const int target_slot = plan.target_slots.front();
  if (target_slot < 0 || target_slot >= static_cast<int>(runtime.slots.size())) {
    return false;
  }
  out_values = runtime.slots[static_cast<std::size_t>(target_slot)];
  return true;
}

bool eval_kernel_query_plan(const KernelProgram &program,
                            const KernelQueryPlan &plan,
                            const KernelEventEvaluator &event_eval,
                            const KernelGuardEvaluator &guard_eval,
                            std::vector<KernelNodeValues> &out_values) {
  thread_local KernelRuntimeState runtime;
  if (!runtime.initialized || runtime.program != &program ||
      runtime.slots.size() != program.ops.size()) {
    reset_kernel_runtime(program, runtime);
  } else {
    invalidate_kernel_runtime_from_slot(runtime, 0);
  }
  bool ok = eval_kernel_query_plan_incremental(program, runtime, plan,
                                               event_eval, guard_eval,
                                               out_values);
  invalidate_kernel_runtime_from_slot(runtime, 0);
  return ok;
}

bool eval_kernel_single_query_plan(const KernelProgram &program,
                                   const KernelQueryPlan &plan,
                                   const KernelEventEvaluator &event_eval,
                                   const KernelGuardEvaluator &guard_eval,
                                   KernelNodeValues &out_values) {
  thread_local KernelRuntimeState runtime;
  if (!runtime.initialized || runtime.program != &program ||
      runtime.slots.size() != program.ops.size()) {
    reset_kernel_runtime(program, runtime);
  } else {
    invalidate_kernel_runtime_from_slot(runtime, 0);
  }
  bool ok = eval_kernel_single_query_plan_incremental(
      program, runtime, plan, event_eval, guard_eval, out_values);
  invalidate_kernel_runtime_from_slot(runtime, 0);
  return ok;
}

TrialParamsSoA build_base_trial_params_soa(const NativeContext &ctx) {
  TrialParamsSoA out;
  out.n_acc = static_cast<int>(ctx.accumulators.size());
  out.dist_code.resize(ctx.accumulators.size(), 0);
  out.onset.resize(ctx.accumulators.size(), 0.0);
  out.q.resize(ctx.accumulators.size(), 0.0);
  out.t0.resize(ctx.accumulators.size(), 0.0);
  out.p1.resize(ctx.accumulators.size(), 0.0);
  out.p2.resize(ctx.accumulators.size(), 0.0);
  out.p3.resize(ctx.accumulators.size(), 0.0);
  out.p4.resize(ctx.accumulators.size(), 0.0);
  out.p5.resize(ctx.accumulators.size(), 0.0);
  out.p6.resize(ctx.accumulators.size(), 0.0);
  out.p7.resize(ctx.accumulators.size(), 0.0);
  out.p8.resize(ctx.accumulators.size(), 0.0);
  for (std::size_t i = 0; i < ctx.accumulators.size(); ++i) {
    const auto &acc = ctx.accumulators[i];
    out.dist_code[i] = acc.dist_cfg.code;
    out.onset[i] = acc.onset;
    out.q[i] = acc.q;
    out.t0[i] = acc.dist_cfg.t0;
    out.p1[i] = acc.dist_cfg.p1;
    out.p2[i] = acc.dist_cfg.p2;
    out.p3[i] = acc.dist_cfg.p3;
    out.p4[i] = acc.dist_cfg.p4;
    out.p5[i] = acc.dist_cfg.p5;
    out.p6[i] = acc.dist_cfg.p6;
    out.p7[i] = acc.dist_cfg.p7;
    out.p8[i] = acc.dist_cfg.p8;
  }
  out.valid = true;
  return out;
}

} // namespace uuber
