#include "kernel_jit.h"

namespace uuber {

KernelJitExecutable compile_kernel_jit(const KernelProgram &program) {
  KernelJitExecutable exec;
  if (!program.valid) {
    return exec;
  }
  exec.valid = true;
  exec.evaluate_node =
      [program](int target_node_idx, const KernelEvalNeed &need,
                const KernelEventEvalFn &event_eval,
                const KernelGuardEvalFn &guard_eval,
                KernelNodeValues &out_values) -> bool {
    return eval_kernel_node(program, target_node_idx, need, event_eval,
                            guard_eval, out_values);
  };
  return exec;
}

} // namespace uuber

