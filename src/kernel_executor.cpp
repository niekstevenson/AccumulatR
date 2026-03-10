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

inline void assign_kernel_batch_defaults(KernelNodeBatchValues &values,
                                         std::size_t point_count) {
  values.density.assign(point_count, 0.0);
  values.survival.assign(point_count, 1.0);
  values.cdf.assign(point_count, 0.0);
}

inline bool kernel_batch_values_match_size(const KernelNodeBatchValues &values,
                                           std::size_t point_count) {
  return values.density.size() == point_count &&
         values.survival.size() == point_count &&
         values.cdf.size() == point_count;
}

inline void clamp_kernel_batch_values(KernelNodeBatchValues &values) {
  const std::size_t point_count = values.density.size();
  for (std::size_t i = 0; i < point_count; ++i) {
    values.density[i] = safe_density(values.density[i]);
    values.survival[i] = clamp_probability(values.survival[i]);
    values.cdf[i] = clamp_probability(values.cdf[i]);
  }
}

inline bool resolve_batch_target_slots(const KernelProgram &program,
                                       const std::vector<int> &target_node_indices,
                                       std::vector<int> &target_slots,
                                       int &max_node_idx) {
  target_slots.clear();
  if (!program.valid || target_node_indices.empty()) {
    return false;
  }
  target_slots.reserve(target_node_indices.size());
  int max_slot = -1;
  for (int node_idx : target_node_indices) {
    if (node_idx < 0 ||
        node_idx >= static_cast<int>(program.outputs.node_idx_to_slot.size())) {
      return false;
    }
    const int slot =
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
  max_node_idx =
      program.outputs.slot_to_node_idx[static_cast<std::size_t>(max_slot)];
  return max_node_idx >= 0;
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

void reset_kernel_batch_runtime(const KernelProgram &program,
                                KernelBatchRuntimeState &runtime,
                                std::size_t point_count) {
  runtime.program = &program;
  runtime.slots.assign(program.ops.size(), KernelNodeBatchValues{});
  const int max_children = std::max(0, program.max_child_count);
  runtime.child_primary.assign(static_cast<std::size_t>(max_children), 0.0);
  runtime.child_density.assign(static_cast<std::size_t>(max_children), 0.0);
  runtime.prefix.assign(static_cast<std::size_t>(max_children + 1), 1.0);
  runtime.suffix.assign(static_cast<std::size_t>(max_children + 1), 1.0);
  runtime.point_count = point_count;
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

void invalidate_kernel_batch_runtime_from_slot(KernelBatchRuntimeState &runtime,
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

bool eval_kernel_node_batch_incremental(
    const KernelProgram &program, KernelBatchRuntimeState &runtime,
    int target_node_idx, const std::vector<double> &times,
    const KernelEvalNeed &need, const KernelEventBatchEvalFn &event_eval_batch,
    const KernelGuardBatchEvalFn &guard_eval_batch,
    KernelNodeBatchValues &out_values) {
  (void)need;
  out_values = KernelNodeBatchValues{};
  if (!program.valid || target_node_idx < 0 || !event_eval_batch) {
    return false;
  }
  if (target_node_idx >=
      static_cast<int>(program.outputs.node_idx_to_slot.size())) {
    return false;
  }
  const int target_slot =
      program.outputs.node_idx_to_slot[static_cast<std::size_t>(target_node_idx)];
  if (target_slot < 0 || target_slot >= static_cast<int>(program.ops.size())) {
    return false;
  }

  const std::size_t point_count = times.size();
  if (!runtime.initialized || runtime.program != &program ||
      runtime.slots.size() != program.ops.size() ||
      runtime.point_count != point_count) {
    reset_kernel_batch_runtime(program, runtime, point_count);
  }
  if (target_slot <= runtime.computed_upto) {
    out_values = runtime.slots[static_cast<std::size_t>(target_slot)];
    return kernel_batch_values_match_size(out_values, point_count);
  }

  const KernelEvalNeed internal_need{true, true, true};
  const int eval_slot_count = target_slot + 1;
  for (int op_idx = std::max(0, runtime.computed_upto + 1); op_idx <= target_slot;
       ++op_idx) {
    const KernelOp &op = program.ops[static_cast<std::size_t>(op_idx)];
    KernelNodeBatchValues out;
    assign_kernel_batch_defaults(out, point_count);
    switch (op.code) {
    case KernelOpCode::Event: {
      if (!event_eval_batch(op.event_idx, times, internal_need, out) ||
          !kernel_batch_values_match_size(out, point_count)) {
        return false;
      }
      clamp_kernel_batch_values(out);
      break;
    }
    case KernelOpCode::And: {
      if (op.child_begin < 0 || op.child_count <= 0 ||
          op.child_begin + op.child_count >
              static_cast<int>(program.children.size())) {
        break;
      }
      const int n = op.child_count;
      for (std::size_t point_idx = 0; point_idx < point_count; ++point_idx) {
        for (int i = 0; i < n; ++i) {
          const int slot =
              program.children[static_cast<std::size_t>(op.child_begin + i)];
          if (slot < 0 || slot >= eval_slot_count) {
            runtime.child_primary[static_cast<std::size_t>(i)] = 0.0;
            runtime.child_density[static_cast<std::size_t>(i)] = 0.0;
            continue;
          }
          const KernelNodeBatchValues &child =
              runtime.slots[static_cast<std::size_t>(slot)];
          runtime.child_primary[static_cast<std::size_t>(i)] =
              clamp_probability(child.cdf[point_idx]);
          runtime.child_density[static_cast<std::size_t>(i)] =
              safe_density(child.density[point_idx]);
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
        out.cdf[point_idx] =
            clamp_probability(runtime.prefix[static_cast<std::size_t>(n)]);
        out.survival[point_idx] = clamp_probability(1.0 - out.cdf[point_idx]);
        double density_sum = 0.0;
        for (int i = 0; i < n; ++i) {
          const double others =
              runtime.prefix[static_cast<std::size_t>(i)] *
              runtime.suffix[static_cast<std::size_t>(i + 1)];
          density_sum += runtime.child_density[static_cast<std::size_t>(i)] *
                         clamp_probability(others);
        }
        out.density[point_idx] = safe_density(density_sum);
      }
      break;
    }
    case KernelOpCode::Or: {
      if (op.child_begin < 0 || op.child_count <= 0 ||
          op.child_begin + op.child_count >
              static_cast<int>(program.children.size())) {
        break;
      }
      const int n = op.child_count;
      for (std::size_t point_idx = 0; point_idx < point_count; ++point_idx) {
        for (int i = 0; i < n; ++i) {
          const int slot =
              program.children[static_cast<std::size_t>(op.child_begin + i)];
          if (slot < 0 || slot >= eval_slot_count) {
            runtime.child_primary[static_cast<std::size_t>(i)] = 1.0;
            runtime.child_density[static_cast<std::size_t>(i)] = 0.0;
            continue;
          }
          const KernelNodeBatchValues &child =
              runtime.slots[static_cast<std::size_t>(slot)];
          runtime.child_primary[static_cast<std::size_t>(i)] =
              clamp_probability(child.survival[point_idx]);
          runtime.child_density[static_cast<std::size_t>(i)] =
              safe_density(child.density[point_idx]);
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
        out.survival[point_idx] =
            clamp_probability(runtime.prefix[static_cast<std::size_t>(n)]);
        out.cdf[point_idx] = clamp_probability(1.0 - out.survival[point_idx]);
        double density_sum = 0.0;
        for (int i = 0; i < n; ++i) {
          const double others =
              runtime.prefix[static_cast<std::size_t>(i)] *
              runtime.suffix[static_cast<std::size_t>(i + 1)];
          density_sum += runtime.child_density[static_cast<std::size_t>(i)] *
                         clamp_probability(others);
        }
        out.density[point_idx] = safe_density(density_sum);
      }
      break;
    }
    case KernelOpCode::Not: {
      if (op.child_begin < 0 || op.child_count <= 0 ||
          op.child_begin >= static_cast<int>(program.children.size())) {
        break;
      }
      const int child_slot =
          program.children[static_cast<std::size_t>(op.child_begin)];
      if (child_slot < 0 || child_slot >= eval_slot_count) {
        break;
      }
      const KernelNodeBatchValues &child =
          runtime.slots[static_cast<std::size_t>(child_slot)];
      for (std::size_t point_idx = 0; point_idx < point_count; ++point_idx) {
        out.cdf[point_idx] = clamp_probability(1.0 - child.cdf[point_idx]);
        out.survival[point_idx] = clamp_probability(child.cdf[point_idx]);
        out.density[point_idx] = 0.0;
      }
      break;
    }
    case KernelOpCode::Guard: {
      KernelNodeBatchValues ref_values;
      KernelNodeBatchValues block_values;
      assign_kernel_batch_defaults(ref_values, point_count);
      assign_kernel_batch_defaults(block_values, point_count);
      if (op.ref_slot >= 0 && op.ref_slot < eval_slot_count) {
        ref_values = runtime.slots[static_cast<std::size_t>(op.ref_slot)];
      }
      if (op.blocker_slot >= 0 && op.blocker_slot < eval_slot_count) {
        block_values =
            runtime.slots[static_cast<std::size_t>(op.blocker_slot)];
      }
      if (guard_eval_batch) {
        if (!guard_eval_batch(op, times, ref_values, block_values, internal_need,
                              out) ||
            !kernel_batch_values_match_size(out, point_count)) {
          return false;
        }
        clamp_kernel_batch_values(out);
      } else {
        for (std::size_t point_idx = 0; point_idx < point_count; ++point_idx) {
          KernelNodeValues ref_value;
          ref_value.density = ref_values.density[point_idx];
          ref_value.survival = ref_values.survival[point_idx];
          ref_value.cdf = ref_values.cdf[point_idx];
          KernelNodeValues block_value;
          block_value.density = block_values.density[point_idx];
          block_value.survival = block_values.survival[point_idx];
          block_value.cdf = block_values.cdf[point_idx];
          const KernelNodeValues point_out =
              default_guard_eval(ref_value, block_value, internal_need);
          out.density[point_idx] = safe_density(point_out.density);
          out.survival[point_idx] = clamp_probability(point_out.survival);
          out.cdf[point_idx] = clamp_probability(point_out.cdf);
        }
      }
      break;
    }
    }
    if (op.out_slot >= 0 &&
        op.out_slot < static_cast<int>(runtime.slots.size())) {
      runtime.slots[static_cast<std::size_t>(op.out_slot)] = std::move(out);
    }
  }

  runtime.computed_upto = std::max(runtime.computed_upto, target_slot);
  out_values = runtime.slots[static_cast<std::size_t>(target_slot)];
  return kernel_batch_values_match_size(out_values, point_count);
}

bool eval_kernel_nodes_batch_incremental(
    const KernelProgram &program, KernelBatchRuntimeState &runtime,
    const std::vector<int> &target_node_indices, const std::vector<double> &times,
    const KernelEvalNeed &need, const KernelEventBatchEvalFn &event_eval_batch,
    const KernelGuardBatchEvalFn &guard_eval_batch,
    std::vector<KernelNodeBatchValues> &out_values) {
  out_values.clear();
  if (!program.valid || target_node_indices.empty() || !event_eval_batch) {
    return false;
  }
  std::vector<int> target_slots;
  int max_node_idx = -1;
  if (!resolve_batch_target_slots(program, target_node_indices, target_slots,
                                  max_node_idx)) {
    return false;
  }
  KernelNodeBatchValues ignored;
  if (!eval_kernel_node_batch_incremental(program, runtime, max_node_idx, times,
                                          need, event_eval_batch,
                                          guard_eval_batch, ignored)) {
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

bool eval_kernel_competitor_product_batch_incremental(
    const KernelProgram &program, KernelBatchRuntimeState &runtime,
    const std::vector<CompetitorCompiledOp> &compiled_ops,
    const std::vector<double> &times,
    const KernelEventBatchEvalFn &event_eval_batch,
    const KernelGuardBatchEvalFn &guard_eval_batch,
    const KernelBatchTransitionApplyFn &apply_transition_batch,
    std::vector<double> &survival_product_out) {
  survival_product_out.assign(times.size(), 1.0);
  if (!program.valid || !event_eval_batch) {
    return false;
  }
  if (compiled_ops.empty()) {
    return true;
  }

  const KernelEvalNeed kernel_need{false, true, true};
  std::vector<int> target_slots;
  for (const CompetitorCompiledOp &op : compiled_ops) {
    if (op.target_node_indices.empty()) {
      if (apply_transition_batch) {
        apply_transition_batch(op, runtime);
      }
      continue;
    }
    int max_node_idx = -1;
    if (!resolve_batch_target_slots(program, op.target_node_indices, target_slots,
                                    max_node_idx)) {
      return false;
    }
    KernelNodeBatchValues ignored;
    if (!eval_kernel_node_batch_incremental(program, runtime, max_node_idx, times,
                                            kernel_need, event_eval_batch,
                                            guard_eval_batch, ignored)) {
      return false;
    }
    for (int slot : target_slots) {
      if (slot < 0 || slot >= static_cast<int>(runtime.slots.size())) {
        return false;
      }
      const KernelNodeBatchValues &vals =
          runtime.slots[static_cast<std::size_t>(slot)];
      if (vals.survival.size() != times.size()) {
        return false;
      }
      for (std::size_t i = 0; i < survival_product_out.size(); ++i) {
        const double surv = clamp_probability(vals.survival[i]);
        if (!std::isfinite(surv) || surv <= 0.0) {
          survival_product_out[i] = 0.0;
        } else if (survival_product_out[i] > 0.0) {
          survival_product_out[i] *= surv;
          if (!std::isfinite(survival_product_out[i]) ||
              survival_product_out[i] <= 0.0) {
            survival_product_out[i] = 0.0;
          }
        }
      }
    }
    if (apply_transition_batch) {
      apply_transition_batch(op, runtime);
    }
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
