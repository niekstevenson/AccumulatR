#pragma once

#include <functional>

#include "kernel_executor.h"

namespace uuber {

// JIT wrapper placeholder: keeps a compiled callable interface so the runtime
// can hard-cutover to a single execution surface without changing APIs.
struct KernelJitExecutable {
  bool valid{false};
  std::function<bool(int, const KernelEvalNeed &, const KernelEventEvalFn &,
                     const KernelGuardEvalFn &, KernelNodeValues &)>
      evaluate_node;
};

KernelJitExecutable compile_kernel_jit(const KernelProgram &program);

} // namespace uuber

