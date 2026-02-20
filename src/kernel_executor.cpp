#include "kernel_executor.h"

#include <cmath>
#include <limits>
#include <vector>

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
                                           const KernelEvalNeed &need) {
  KernelNodeValues out;
  out.density = need.density ? safe_density(ref.density * block.survival) : 0.0;
  if (need.cdf) {
    out.cdf = clamp_probability(ref.cdf * block.survival);
  }
  if (need.survival) {
    out.survival = need.cdf ? clamp_probability(1.0 - out.cdf) : block.survival;
  }
  return out;
}

} // namespace

void reset_kernel_runtime(const KernelProgram &program,
                          KernelRuntimeState &runtime) {
  runtime.program = &program;
  runtime.slots.assign(program.ops.size(), KernelNodeValues{});
  const int max_children = std::max(0, program.max_child_count);
  runtime.child_primary.assign(static_cast<std::size_t>(max_children), 0.0);
  runtime.child_density.assign(static_cast<std::size_t>(max_children), 0.0);
  runtime.prefix.assign(static_cast<std::size_t>(max_children + 1), 1.0);
  runtime.suffix.assign(static_cast<std::size_t>(max_children + 1), 1.0);
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

bool eval_kernel_node_incremental(const KernelProgram &program,
                                  KernelRuntimeState &runtime,
                                  int target_node_idx,
                                  const KernelEvalNeed &need,
                                  const KernelEventEvalFn &event_eval,
                                  const KernelGuardEvalFn &guard_eval,
                                  KernelNodeValues &out_values) {
  if (!program.valid || target_node_idx < 0 || !event_eval) {
    return false;
  }
  if (target_node_idx >=
      static_cast<int>(program.outputs.node_idx_to_slot.size())) {
    return false;
  }
  int target_slot =
      program.outputs.node_idx_to_slot[static_cast<std::size_t>(target_node_idx)];
  if (target_slot < 0 || target_slot >= static_cast<int>(program.ops.size())) {
    return false;
  }
  if (!runtime.initialized || runtime.program != &program ||
      runtime.slots.size() != program.ops.size()) {
    reset_kernel_runtime(program, runtime);
  }
  if (target_slot <= runtime.computed_upto) {
    out_values = runtime.slots[static_cast<std::size_t>(target_slot)];
    return true;
  }

  const int eval_slot_count = target_slot + 1;
  for (int op_idx = std::max(0, runtime.computed_upto + 1); op_idx <= target_slot;
       ++op_idx) {
    const KernelOp &op = program.ops[static_cast<std::size_t>(op_idx)];
    KernelNodeValues out;
    switch (op.code) {
    case KernelOpCode::Event: {
      out = event_eval(op.event_idx);
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
            clamp_probability(runtime.slots[slot].cdf);
        runtime.child_density[static_cast<std::size_t>(i)] =
            safe_density(runtime.slots[slot].density);
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
      out.cdf =
          clamp_probability(runtime.prefix[static_cast<std::size_t>(n)]);
      if (need.survival) {
        out.survival = clamp_probability(1.0 - out.cdf);
      }
      if (need.density) {
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
            clamp_probability(runtime.slots[slot].survival);
        runtime.child_density[static_cast<std::size_t>(i)] =
            safe_density(runtime.slots[slot].density);
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
      out.survival =
          clamp_probability(runtime.prefix[static_cast<std::size_t>(n)]);
      out.cdf = clamp_probability(1.0 - out.survival);
      if (need.density) {
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
      out.cdf = clamp_probability(1.0 - child.cdf);
      out.survival = clamp_probability(child.cdf);
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
        out = guard_eval(op, ref, block, need);
      } else {
        out = default_guard_eval(ref, block, need);
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

  runtime.computed_upto = std::max(runtime.computed_upto, target_slot);
  out_values = runtime.slots[static_cast<std::size_t>(target_slot)];
  return true;
}

bool eval_kernel_nodes_incremental(const KernelProgram &program,
                                   KernelRuntimeState &runtime,
                                   const std::vector<int> &target_node_indices,
                                   const KernelEvalNeed &need,
                                   const KernelEventEvalFn &event_eval,
                                   const KernelGuardEvalFn &guard_eval,
                                   std::vector<KernelNodeValues> &out_values) {
  out_values.clear();
  if (!program.valid || target_node_indices.empty() || !event_eval) {
    return false;
  }
  std::vector<int> target_slots;
  target_slots.reserve(target_node_indices.size());
  int max_slot = -1;
  for (int node_idx : target_node_indices) {
    if (node_idx < 0 ||
        node_idx >= static_cast<int>(program.outputs.node_idx_to_slot.size())) {
      return false;
    }
    int slot =
        program.outputs.node_idx_to_slot[static_cast<std::size_t>(node_idx)];
    if (slot < 0 || slot >= static_cast<int>(program.ops.size())) {
      return false;
    }
    target_slots.push_back(slot);
    if (slot > max_slot) {
      max_slot = slot;
    }
  }
  if (max_slot < 0 ||
      max_slot >= static_cast<int>(program.outputs.slot_to_node_idx.size())) {
    return false;
  }
  const int max_node_idx =
      program.outputs.slot_to_node_idx[static_cast<std::size_t>(max_slot)];
  KernelNodeValues ignored{};
  if (!eval_kernel_node_incremental(program, runtime, max_node_idx, need,
                                    event_eval, guard_eval, ignored)) {
    return false;
  }
  out_values.reserve(target_slots.size());
  for (int slot : target_slots) {
    if (slot < 0 || slot >= static_cast<int>(runtime.slots.size())) {
      return false;
    }
    out_values.push_back(runtime.slots[static_cast<std::size_t>(slot)]);
  }
  return true;
}

bool eval_kernel_node(const KernelProgram &program, int target_node_idx,
                      const KernelEvalNeed &need,
                      const KernelEventEvalFn &event_eval,
                      const KernelGuardEvalFn &guard_eval,
                      KernelNodeValues &out_values) {
  thread_local KernelRuntimeState runtime;
  if (!runtime.initialized || runtime.program != &program ||
      runtime.slots.size() != program.ops.size()) {
    reset_kernel_runtime(program, runtime);
  } else {
    invalidate_kernel_runtime_from_slot(runtime, 0);
  }
  bool ok = eval_kernel_node_incremental(program, runtime, target_node_idx, need,
                                         event_eval, guard_eval, out_values);
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
  for (std::size_t i = 0; i < ctx.accumulators.size(); ++i) {
    const auto &acc = ctx.accumulators[i];
    out.dist_code[i] = acc.dist_cfg.code;
    out.onset[i] = acc.onset;
    out.q[i] = acc.q;
    out.t0[i] = acc.dist_cfg.t0;
    out.p1[i] = acc.dist_cfg.p1;
    out.p2[i] = acc.dist_cfg.p2;
    out.p3[i] = acc.dist_cfg.p3;
  }
  out.valid = true;
  return out;
}

} // namespace uuber
