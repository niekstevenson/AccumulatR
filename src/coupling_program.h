#pragma once

#include <vector>

#include "vector_runtime.h"

uuber::VectorProgram resolve_outcome_coupling_program_with_generic(
    const uuber::NativeContext &ctx, int node_id,
    const std::vector<int> &competitor_node_ids, bool forced_empty);
