#pragma once

#include "exact_common_types.hpp"

namespace accumulatr::eval {
namespace detail {

struct ExactComplexityMetrics {
  semantic::Index symbolic_region_count{0};
  semantic::Index symbolic_cell_count{0};
  semantic::Index max_symbolic_cells_per_region{0};
  semantic::Index negative_symbolic_cell_count{0};
  semantic::Index overlapping_symbolic_cell_pair_count{0};
  semantic::Index expr_relation_atom_count{0};
  semantic::Index compiled_root_count{0};
  semantic::Index compiled_node_count{0};
  semantic::Index integral_node_count{0};
  semantic::Index integral_kernel_count{0};
  semantic::Index source_product_integral_kernel_count{0};
  semantic::Index generic_integral_kernel_count{0};
  semantic::Index max_integral_depth{0};
};

} // namespace detail
} // namespace accumulatr::eval
