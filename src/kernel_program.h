#pragma once

#include "kernel_executor.h"

namespace uuber {

// Compile dense IR nodes to a single topologically ordered opcode program.
KernelProgram compile_kernel_program(const IrContext &ir);

const KernelQueryPlan *precompute_kernel_query_plan(KernelProgram &program,
                                                    int target_node_idx,
                                                    EvalNeed need);

int precompute_kernel_multi_query_plan(KernelProgram &program,
                                       const std::vector<int> &target_node_indices,
                                       const std::vector<EvalNeed> &target_needs);

const KernelQueryPlan *get_kernel_query_plan(const KernelProgram &program,
                                             int target_node_idx,
                                             EvalNeed need);

const KernelQueryPlan *get_kernel_multi_query_plan(const KernelProgram &program,
                                                   int query_plan_index);

} // namespace uuber
