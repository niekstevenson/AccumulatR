#pragma once

#include <vector>

#include "context.h"

const uuber::CompetitorClusterCacheEntry &
fetch_competitor_cluster_cache(const uuber::NativeContext &ctx,
                               const std::vector<int> &competitor_ids);
