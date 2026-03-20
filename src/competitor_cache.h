#pragma once

#include <vector>

#include "context.h"

struct NodeEvalState;

namespace uuber {
struct KernelBatchRuntimeState;
}

const uuber::CompetitorClusterCacheEntry &
fetch_competitor_cluster_cache(const uuber::NativeContext &ctx,
                               const std::vector<int> &competitor_ids);

void competitor_survival_batch_from_state_compiled_ops(
    const uuber::NativeContext &ctx, const std::vector<int> &competitor_ids,
    NodeEvalState &state, const std::vector<double> &times,
    std::vector<double> &survival_out,
    uuber::KernelBatchRuntimeState *kernel_batch_runtime = nullptr,
    const uuber::CompetitorClusterCacheEntry *competitor_cache = nullptr);
