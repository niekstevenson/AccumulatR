#pragma once

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

#include "compiled_math_types.hpp"
#include "exact_common_types.hpp"
#include "../runtime/exact_program.hpp"

namespace accumulatr::eval {
namespace detail {

struct ExactTriggerState {
  double weight{1.0};
  const std::uint8_t *shared_started{nullptr};
};

struct ExactCompiledTriggerWeightTerm {
  semantic::Index leaf_index{semantic::kInvalidIndex};
  std::uint8_t shared_started{2U};
};

struct ExactCompiledTriggerState {
  double fixed_weight{1.0};
  ExactIndexSpan weight_terms{};
  semantic::Index shared_started_offset{0};
};

struct ExactCompiledTriggerStateTable {
  std::vector<ExactCompiledTriggerState> states;
  std::vector<ExactCompiledTriggerWeightTerm> weight_terms;
  std::vector<std::uint8_t> shared_started_values;
  semantic::Index trigger_count{0};
};

struct ExactSequenceState {
  double lower_bound{0.0};
  std::vector<double> exact_times;
  std::vector<double> upper_bounds;
  std::vector<double> expr_upper_bounds;
  std::vector<double> expr_upper_normalizers;
};

struct ExactRankedFrontierEntry {
  double probability{0.0};
  semantic::Index state_index{semantic::kInvalidIndex};
};

struct ExactStepDistributionView {
  double total_probability{0.0};
  const std::vector<double> *transition_probabilities{nullptr};
};

struct ExactExprKernel {
  semantic::ExprKind kind{semantic::ExprKind::Impossible};
  ExactIndexSpan children;
  semantic::Index event_source_id{semantic::kInvalidIndex};
  semantic::Index guard_ref_expr_id{semantic::kInvalidIndex};
  semantic::Index guard_blocker_expr_id{semantic::kInvalidIndex};
  bool has_unless{false};
};

struct ExactSourceKernel {
  CompiledSourceChannelKernelKind kind{
      CompiledSourceChannelKernelKind::Invalid};
  semantic::Index source_id{semantic::kInvalidIndex};
  semantic::Index leaf_index{semantic::kInvalidIndex};
  semantic::Index pool_index{semantic::kInvalidIndex};
  semantic::Index onset_source_id{semantic::kInvalidIndex};
  semantic::Index pool_member_offset{0};
  semantic::Index pool_member_count{0};
  semantic::Index pool_k{0};
};

inline std::vector<semantic::Index> merge_sorted_support(
    std::vector<semantic::Index> merged,
    const std::vector<semantic::Index> &rhs) {
  for (const auto idx : rhs) {
    if (std::find(merged.begin(), merged.end(), idx) == merged.end()) {
      merged.push_back(idx);
    }
  }
  std::sort(merged.begin(), merged.end());
  return merged;
}

inline bool supports_overlap(const std::vector<semantic::Index> &lhs,
                             const std::vector<semantic::Index> &rhs) {
  std::vector<semantic::Index> overlap;
  std::set_intersection(lhs.begin(),
                        lhs.end(),
                        rhs.begin(),
                        rhs.end(),
                        std::back_inserter(overlap));
  return !overlap.empty();
}

inline bool support_contains_source(const std::vector<semantic::Index> &support,
                                    const semantic::Index source_id) {
  return source_id != semantic::kInvalidIndex &&
         std::binary_search(support.begin(), support.end(), source_id);
}


inline bool has_reason(const std::vector<std::string> &reasons,
                       const std::string &needle) {
  return std::find(reasons.begin(), reasons.end(), needle) != reasons.end();
}

class ExactSupportBuilder {
public:
  explicit ExactSupportBuilder(const runtime::LoweredExactVariant &lowered)
      : program_(lowered.program),
        leaf_ready_(static_cast<std::size_t>(program_.layout.n_leaves), 0U),
        pool_ready_(static_cast<std::size_t>(program_.layout.n_pools), 0U),
        expr_ready_(program_.expr_kind.size(), 0U),
        leaf_supports_(static_cast<std::size_t>(program_.layout.n_leaves)),
        pool_supports_(static_cast<std::size_t>(program_.layout.n_pools)),
        expr_supports_(program_.expr_kind.size()) {}

  std::vector<std::vector<semantic::Index>> build_leaf_supports() {
    for (semantic::Index i = 0; i < program_.layout.n_leaves; ++i) {
      leaf_support(i);
    }
    return leaf_supports_;
  }

  std::vector<std::vector<semantic::Index>> build_pool_supports() {
    for (semantic::Index i = 0; i < program_.layout.n_pools; ++i) {
      pool_support(i);
    }
    return pool_supports_;
  }

  std::vector<std::vector<semantic::Index>> build_expr_supports() {
    for (semantic::Index i = 0;
         i < static_cast<semantic::Index>(program_.expr_kind.size());
         ++i) {
      expr_support(i);
    }
    return expr_supports_;
  }

private:
  const runtime::ExactProgram &program_;
  std::vector<std::uint8_t> leaf_ready_;
  std::vector<std::uint8_t> pool_ready_;
  std::vector<std::uint8_t> expr_ready_;
  std::vector<std::vector<semantic::Index>> leaf_supports_;
  std::vector<std::vector<semantic::Index>> pool_supports_;
  std::vector<std::vector<semantic::Index>> expr_supports_;

  const std::vector<semantic::Index> &leaf_support(const semantic::Index idx) {
    const auto pos = static_cast<std::size_t>(idx);
    if (leaf_ready_[pos] == 2U) {
      return leaf_supports_[pos];
    }
    if (leaf_ready_[pos] == 1U) {
      throw std::runtime_error("cyclic leaf dependency in exact support builder");
    }
    leaf_ready_[pos] = 1U;
    std::vector<semantic::Index> support{idx};
    const auto onset_kind =
        static_cast<semantic::OnsetKind>(program_.onset_kind[pos]);
    if (onset_kind != semantic::OnsetKind::Absolute) {
      support = merge_sorted_support(
          std::move(support),
          source_support(
              static_cast<semantic::SourceKind>(program_.onset_source_kind[pos]),
              program_.onset_source_index[pos]));
    }
    leaf_supports_[pos] = std::move(support);
    leaf_ready_[pos] = 2U;
    return leaf_supports_[pos];
  }

  const std::vector<semantic::Index> &pool_support(const semantic::Index idx) {
    const auto pos = static_cast<std::size_t>(idx);
    if (pool_ready_[pos] == 2U) {
      return pool_supports_[pos];
    }
    if (pool_ready_[pos] == 1U) {
      throw std::runtime_error("cyclic pool dependency in exact support builder");
    }
    pool_ready_[pos] = 1U;
    std::vector<semantic::Index> support;
    const auto begin = program_.pool_member_offsets[pos];
    const auto end = program_.pool_member_offsets[pos + 1U];
    for (semantic::Index i = begin; i < end; ++i) {
      support = merge_sorted_support(
          std::move(support),
          source_support(
              static_cast<semantic::SourceKind>(
                  program_.pool_member_kind[static_cast<std::size_t>(i)]),
              program_.pool_member_indices[static_cast<std::size_t>(i)]));
    }
    pool_supports_[pos] = std::move(support);
    pool_ready_[pos] = 2U;
    return pool_supports_[pos];
  }

  const std::vector<semantic::Index> &expr_support(const semantic::Index idx) {
    const auto pos = static_cast<std::size_t>(idx);
    if (expr_ready_[pos] == 2U) {
      return expr_supports_[pos];
    }
    if (expr_ready_[pos] == 1U) {
      throw std::runtime_error("cyclic expression dependency in exact support builder");
    }
    expr_ready_[pos] = 1U;
    std::vector<semantic::Index> support;
    const auto kind =
        static_cast<semantic::ExprKind>(program_.expr_kind[pos]);
    switch (kind) {
    case semantic::ExprKind::Event:
      support = source_support(
          static_cast<semantic::SourceKind>(program_.expr_source_kind[pos]),
          program_.expr_source_index[pos]);
      break;
    case semantic::ExprKind::And:
    case semantic::ExprKind::Or:
    case semantic::ExprKind::Not: {
      const auto begin = program_.expr_arg_offsets[pos];
      const auto end = program_.expr_arg_offsets[pos + 1U];
      for (semantic::Index i = begin; i < end; ++i) {
        support = merge_sorted_support(
            std::move(support),
            expr_support(program_.expr_args[static_cast<std::size_t>(i)]));
      }
      break;
    }
    case semantic::ExprKind::Guard:
      support = merge_sorted_support(
          std::move(support),
          expr_support(program_.expr_ref_child[pos]));
      support = merge_sorted_support(
          std::move(support),
          expr_support(program_.expr_blocker_child[pos]));
      for (semantic::Index i = program_.expr_arg_offsets[pos];
           i < program_.expr_arg_offsets[pos + 1U];
           ++i) {
        support = merge_sorted_support(
            std::move(support),
            expr_support(program_.expr_args[static_cast<std::size_t>(i)]));
      }
      break;
    case semantic::ExprKind::Impossible:
    case semantic::ExprKind::TrueExpr:
      break;
    }
    expr_supports_[pos] = std::move(support);
    expr_ready_[pos] = 2U;
    return expr_supports_[pos];
  }

  std::vector<semantic::Index> source_support(const semantic::SourceKind kind,
                                              const semantic::Index index) {
    switch (kind) {
    case semantic::SourceKind::Leaf:
      return leaf_support(index);
    case semantic::SourceKind::Pool:
      return pool_support(index);
    case semantic::SourceKind::Special:
      break;
    }
    return {};
  }
};

inline semantic::Index child_event_source_index(const runtime::ExactProgram &program,
                                                const semantic::Index expr_idx) {
  return program.expr_source_index[static_cast<std::size_t>(expr_idx)];
}

inline semantic::SourceKind child_event_source_kind(
    const runtime::ExactProgram &program,
    const semantic::Index expr_idx) {
  return static_cast<semantic::SourceKind>(
      program.expr_source_kind[static_cast<std::size_t>(expr_idx)]);
}

inline void validate_exact_expr(const runtime::LoweredExactVariant &lowered,
                                const semantic::Index expr_idx) {
  const auto &program = lowered.program;
  const auto kind = static_cast<semantic::ExprKind>(
      program.expr_kind[static_cast<std::size_t>(expr_idx)]);
  if (kind == semantic::ExprKind::Event) {
    return;
  }
  if (kind == semantic::ExprKind::Impossible || kind == semantic::ExprKind::TrueExpr) {
    return;
  }
  if (kind == semantic::ExprKind::Not) {
    validate_exact_expr(
        lowered,
        program.expr_args[static_cast<std::size_t>(
            program.expr_arg_offsets[static_cast<std::size_t>(expr_idx)])]);
    return;
  }
  if (kind == semantic::ExprKind::And || kind == semantic::ExprKind::Or) {
    for (semantic::Index i = program.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
         i < program.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
         ++i) {
      validate_exact_expr(lowered, program.expr_args[static_cast<std::size_t>(i)]);
    }
    return;
  }
  if (kind == semantic::ExprKind::Guard) {
    validate_exact_expr(
        lowered,
        program.expr_ref_child[static_cast<std::size_t>(expr_idx)]);
    validate_exact_expr(
        lowered,
        program.expr_blocker_child[static_cast<std::size_t>(expr_idx)]);
    for (semantic::Index i = program.expr_arg_offsets[static_cast<std::size_t>(expr_idx)];
         i < program.expr_arg_offsets[static_cast<std::size_t>(expr_idx + 1)];
         ++i) {
      validate_exact_expr(lowered, program.expr_args[static_cast<std::size_t>(i)]);
    }
  }
}

} // namespace detail
} // namespace accumulatr::eval
