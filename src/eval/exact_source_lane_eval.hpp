#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

#include "compiled_lane_workspace.hpp"
#include "exact_source_leaf_batch.hpp"

namespace accumulatr::eval {
namespace detail {

struct SourceLaneScratch {
  void ensure_items(const std::size_t count) {
    lanes.resize(count);
    lanes_b.resize(count);
    lanes_c.resize(count);
    positions.resize(count);
    positions_b.resize(count);
    positions_c.resize(count);
    times.resize(count);
    times_b.resize(count);
    times_c.resize(count);
    parent_times.resize(count);
    weights.resize(count);
    lower_bounds.resize(count);
    upper_bounds.resize(count);
    exact_flags.resize(count);
  }

  void ensure_pool(const std::size_t count,
                   const std::size_t members,
                   const bool need_pdf) {
    weights.resize(count);
    const auto member_values = count * members;
    if (need_pdf) {
      member_pdf.resize(member_values);
    }
    member_cdf.resize(member_values);
    member_survival.resize(member_values);
    const auto width = members + 1U;
    const auto table_values = count * width * width;
    prefix.resize(table_values);
    if (need_pdf) {
      suffix.resize(table_values);
    }
  }

  std::vector<semantic::Index> lanes;
  std::vector<semantic::Index> lanes_b;
  std::vector<semantic::Index> lanes_c;
  std::vector<std::size_t> positions;
  std::vector<std::size_t> positions_b;
  std::vector<std::size_t> positions_c;
  std::vector<double> times;
  std::vector<double> times_b;
  std::vector<double> times_c;
  std::vector<double> parent_times;
  std::vector<double> weights;
  std::vector<double> lower_bounds;
  std::vector<double> upper_bounds;
  std::vector<std::uint8_t> exact_flags;
  std::vector<double> cdf_totals;
  PreparedSourceLeafBatch leaf_batch;
  SourceLaneFill out;
  SourceLaneFill unconditioned;
  SourceLaneFill lower;
  SourceLaneFill upper;
  SourceLaneFill shifted;
  std::vector<double> member_pdf;
  std::vector<double> member_cdf;
  std::vector<double> member_survival;
  std::vector<double> prefix;
  std::vector<double> suffix;
};

struct SourceLaneWorkspace {
  SourceLaneScratch &layer(const std::size_t depth) {
    while (layers.size() <= depth) {
      layers.push_back(std::make_unique<SourceLaneScratch>());
    }
    return *layers[depth];
  }

  std::vector<std::unique_ptr<SourceLaneScratch>> layers;
};

inline void source_lane_fill_forced(
    SourceLaneFill *out,
    const std::size_t position,
    const ExactRelation relation) {
  source_lane_store_fill(
      out,
      position,
      exact_source_forced_fill(relation, out->mask));
}

template <leaf::DistKind Kind,
          bool IdentitySourceLanes,
          typename ElapsedTime>
inline void prepare_source_leaf_batch(
    const CompiledMathSourceProductProgram &source_program,
    CompiledLaneWorkspace *workspace,
    const CompiledLaneFrame &frame,
    const semantic::Index *lanes,
    const std::size_t count,
    ElapsedTime &elapsed_time,
    PreparedSourceLeafBatch *out) {
  constexpr auto parameter_count = static_cast<std::size_t>(
      leaf::dist_param_count(Kind));
  out->ensure_lanes(count);
  const auto leaves = workspace->source_state().leaf_batch(
      source_program.leaf_index);
  auto *elapsed = out->elapsed();
  auto *q = out->q();
  std::array<double *, leaf::kMaxDistParamCount> parameters{};
  for (std::size_t slot = 0U; slot < parameter_count; ++slot) {
    parameters[slot] = out->parameter(slot);
  }
  for (std::size_t i = 0U; i < count; ++i) {
    std::size_t source_lane = i;
    if constexpr (!IdentitySourceLanes) {
      const auto frame_lane = static_cast<std::size_t>(
          compiled_frame_lane(lanes, i));
      source_lane = workspace->source_lane(frame, frame_lane);
    }
    const int physical_row = leaves.physical_row(source_lane);
    elapsed[i] = elapsed_time(i, leaves, physical_row);
    q[i] = leaves.q(source_lane, physical_row);
    for (std::size_t slot = 0U; slot < parameter_count; ++slot) {
      parameters[slot][i] = leaves.param(slot, physical_row);
    }
  }
}

template <leaf::DistKind Kind, typename ElapsedTime>
inline void source_lane_leaf_program_kind(
    const CompiledMathSourceProductProgram &source_program,
    CompiledLaneWorkspace *workspace,
    const CompiledLaneFrame &frame,
    const semantic::Index *lanes,
    const std::size_t count,
    ElapsedTime &elapsed_time,
    PreparedSourceLeafBatch *leaf_batch,
    SourceLaneFill *out) {
  if (lanes == nullptr && frame.source_lanes_identity) {
    prepare_source_leaf_batch<Kind, true>(
        source_program, workspace, frame, lanes, count,
        elapsed_time, leaf_batch);
  } else {
    prepare_source_leaf_batch<Kind, false>(
        source_program, workspace, frame, lanes, count,
        elapsed_time, leaf_batch);
  }
  evaluate_prepared_source_leaf_batch<Kind>(*leaf_batch, out);
}

template <typename ElapsedTime>
inline void source_lane_leaf_program(
    const CompiledMathSourceProductProgram &source_program,
    CompiledLaneWorkspace *workspace,
    const CompiledLaneFrame &frame,
    const semantic::Index *lanes,
    const std::size_t count,
    ElapsedTime &&elapsed_time,
    PreparedSourceLeafBatch *leaf_batch,
    SourceLaneFill *out) {
  switch (static_cast<leaf::DistKind>(source_program.leaf_dist_kind)) {
  case leaf::DistKind::Lognormal:
    source_lane_leaf_program_kind<leaf::DistKind::Lognormal>(
        source_program, workspace, frame, lanes, count,
        elapsed_time, leaf_batch, out);
    break;
  case leaf::DistKind::Gamma:
    source_lane_leaf_program_kind<leaf::DistKind::Gamma>(
        source_program, workspace, frame, lanes, count,
        elapsed_time, leaf_batch, out);
    break;
  case leaf::DistKind::Exgauss:
    source_lane_leaf_program_kind<leaf::DistKind::Exgauss>(
        source_program, workspace, frame, lanes, count,
        elapsed_time, leaf_batch, out);
    break;
  case leaf::DistKind::LBA:
    source_lane_leaf_program_kind<leaf::DistKind::LBA>(
        source_program, workspace, frame, lanes, count,
        elapsed_time, leaf_batch, out);
    break;
  case leaf::DistKind::RDM:
    source_lane_leaf_program_kind<leaf::DistKind::RDM>(
        source_program, workspace, frame, lanes, count,
        elapsed_time, leaf_batch, out);
    break;
  }
}

inline SourceLaneFill &evaluate_source_program_lanes(
    const CompiledMathProgram &program,
    const semantic::Index source_program_id,
    CompiledLaneWorkspace *workspace,
    const CompiledLaneFrame &frame,
    const semantic::Index *lanes,
    const double *times,
    const std::size_t count,
    const std::uint8_t fill_mask,
    SourceLaneWorkspace *source_workspace,
    const std::size_t depth);

inline SourceLaneFill &evaluate_conditioned_source_program_lanes(
    const CompiledMathProgram &program,
    const CompiledMathSourceProductProgram &source_program,
    CompiledLaneWorkspace *workspace,
    const CompiledLaneFrame &frame,
    const semantic::Index *lanes,
    const double *times,
    const std::size_t count,
    const std::uint8_t fill_mask,
    SourceLaneWorkspace *source_workspace,
    const std::size_t depth) {
  auto &scratch = source_workspace->layer(depth);
  scratch.out.assign(count, fill_mask);
  scratch.ensure_items(count);
  std::size_t unresolved = 0U;
  bool any_finite_upper = false;
  for (std::size_t i = 0; i < count; ++i) {
    const auto lane = static_cast<std::size_t>(
        compiled_frame_lane(lanes, i));
    const auto source_lane = workspace->source_lane(frame, lane);
    const double exact = workspace->source_state().sequence_exact_time(
        source_lane, source_program.source_id);
    const double upper = workspace->source_state().sequence_upper_bound(
        source_lane, source_program.source_id);
    const double lower = std::isfinite(exact) || std::isfinite(upper)
                             ? 0.0
                             : workspace->source_state().sequence_lower_bound(
                                   source_lane);
    scratch.lower_bounds[i] = lower;
    scratch.upper_bounds[i] = upper;
    if (std::isfinite(upper) && !(upper > lower)) {
      continue;
    }
    if (std::isfinite(exact)) {
      const bool valid =
          (!(lower > 0.0) || exact > lower) &&
          (!std::isfinite(upper) || exact < upper) &&
          times[i] >= exact;
      if (valid) {
        source_lane_store_fill(
            &scratch.out,
            i,
            exact_source_certain_fill(fill_mask));
      }
      continue;
    }
    scratch.lanes[unresolved] = compiled_frame_lane(lanes, i);
    scratch.positions[unresolved] = i;
    scratch.times[unresolved] = times[i];
    any_finite_upper = any_finite_upper || std::isfinite(upper);
    ++unresolved;
  }
  if (unresolved == 0U) {
    return scratch.out;
  }

  const std::uint8_t unconditioned_mask =
      fill_mask |
      (any_finite_upper && (fill_mask & kLeafChannelSurvival) != 0U
           ? kLeafChannelCdf
           : 0U);
  const auto &unconditioned = evaluate_source_program_lanes(
      program,
      source_program.child_program_id,
      workspace,
      frame,
      scratch.lanes.data(),
      scratch.times.data(),
      unresolved,
      unconditioned_mask,
      source_workspace,
      depth + 1U);
  scratch.unconditioned.resize(unresolved, unconditioned_mask);
  if ((unconditioned_mask & kLeafChannelPdf) != 0U) {
    scratch.unconditioned.pdf = unconditioned.pdf;
  }
  if ((unconditioned_mask & kLeafChannelCdf) != 0U) {
    scratch.unconditioned.cdf = unconditioned.cdf;
  }
  if ((unconditioned_mask & kLeafChannelSurvival) != 0U) {
    scratch.unconditioned.survival = unconditioned.survival;
  }

  std::size_t lower_count = 0U;
  for (std::size_t j = 0; j < unresolved; ++j) {
    const auto position = scratch.positions[j];
    if (scratch.lower_bounds[position] > 0.0) {
      scratch.lanes_b[lower_count] = compiled_frame_lane(lanes, position);
      scratch.positions_b[lower_count] = j;
      scratch.times_b[lower_count] = scratch.lower_bounds[position];
      ++lower_count;
    }
  }
  scratch.lower.resize(unresolved, kLeafChannelCdf | kLeafChannelSurvival);
  if (lower_count > 0U) {
    const auto &lower = evaluate_source_program_lanes(
        program,
        source_program.child_program_id,
        workspace,
        frame,
        scratch.lanes_b.data(),
        scratch.times_b.data(),
        lower_count,
        kLeafChannelCdf | kLeafChannelSurvival,
        source_workspace,
        depth + 1U);
    for (std::size_t k = 0; k < lower_count; ++k) {
      const auto j = scratch.positions_b[k];
      scratch.lower.cdf[j] = lower.cdf[k];
      scratch.lower.survival[j] = lower.survival[k];
    }
  }

  std::size_t upper_count = 0U;
  for (std::size_t j = 0; j < unresolved; ++j) {
    const auto position = scratch.positions[j];
    if (std::isfinite(scratch.upper_bounds[position])) {
      scratch.lanes_c[upper_count] = compiled_frame_lane(lanes, position);
      scratch.positions_c[upper_count] = j;
      scratch.times_c[upper_count] = scratch.upper_bounds[position];
      ++upper_count;
    }
  }
  scratch.upper.resize(unresolved, kLeafChannelCdf);
  if (upper_count > 0U) {
    const auto &upper = evaluate_source_program_lanes(
        program,
        source_program.child_program_id,
        workspace,
        frame,
        scratch.lanes_c.data(),
        scratch.times_c.data(),
        upper_count,
        kLeafChannelCdf,
        source_workspace,
        depth + 1U);
    for (std::size_t k = 0; k < upper_count; ++k) {
      scratch.upper.cdf[scratch.positions_c[k]] = upper.cdf[k];
    }
  }

  for (std::size_t j = 0; j < unresolved; ++j) {
    const auto position = scratch.positions[j];
    const double lower_bound = scratch.lower_bounds[position];
    const double upper_bound = scratch.upper_bounds[position];
    ExactSourceFill raw;
    if ((unconditioned_mask & kLeafChannelPdf) != 0U) {
      raw.pdf = scratch.unconditioned.pdf[j];
    }
    if ((unconditioned_mask & kLeafChannelCdf) != 0U) {
      raw.cdf = scratch.unconditioned.cdf[j];
    }
    if ((unconditioned_mask & kLeafChannelSurvival) != 0U) {
      raw.survival = scratch.unconditioned.survival[j];
    }
    if (!(lower_bound > 0.0) && !std::isfinite(upper_bound)) {
      source_lane_store_fill(&scratch.out, position, raw);
      continue;
    }
    ExactSourceFill lower;
    if (lower_bound > 0.0) {
      lower.cdf = scratch.lower.cdf[j];
      lower.survival = scratch.lower.survival[j];
    }
    if (!std::isfinite(upper_bound)) {
      source_lane_store_fill(
          &scratch.out,
          position,
          exact_source_conditionalize(
              raw, lower, fill_mask));
      continue;
    }
    ExactSourceFill upper;
    upper.cdf = scratch.upper.cdf[j];
    if (!std::isfinite(upper.cdf - lower.cdf) ||
        !(upper.cdf - lower.cdf > 0.0)) {
      continue;
    }
    if (times[position] >= upper_bound) {
      source_lane_store_fill(
          &scratch.out,
          position,
          exact_source_certain_fill(fill_mask));
    } else if (times[position] > lower_bound) {
      source_lane_store_fill(
          &scratch.out,
          position,
          exact_source_conditionalize_between(
              raw, lower, upper, fill_mask));
    }
  }
  return scratch.out;
}

inline SourceLaneFill &evaluate_onset_source_program_lanes(
    const CompiledMathProgram &program,
    const CompiledMathSourceProductProgram &source_program,
    CompiledLaneWorkspace *workspace,
    const CompiledLaneFrame &frame,
    const semantic::Index *lanes,
    const double *times,
    const std::size_t count,
    const std::uint8_t fill_mask,
    SourceLaneWorkspace *source_workspace,
    const std::size_t depth) {
  auto &scratch = source_workspace->layer(depth);
  scratch.out.assign(count, fill_mask);
  const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
  const bool need_cdf =
      (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  const auto shifted_mask = static_cast<std::uint8_t>(
      (need_pdf ? kLeafChannelPdf : 0U) |
      (need_cdf ? kLeafChannelCdf : 0U));

  scratch.ensure_items(kExactExpandedLaneTileSize);
  if (need_cdf) {
    scratch.cdf_totals.resize(count);
    std::fill(
        scratch.cdf_totals.begin(),
        scratch.cdf_totals.end(),
        std::numeric_limits<double>::quiet_NaN());
  }
  constexpr auto parent_tile =
      exact_expanded_parent_lane_tile_size<
          quadrature::kDefaultFiniteOrder>();
  for (std::size_t parent_begin = 0U;
       parent_begin < count;
       parent_begin += parent_tile) {
    const auto parent_end = std::min(count, parent_begin + parent_tile);
    std::size_t child_count = 0U;
    for (std::size_t i = parent_begin; i < parent_end; ++i) {
      const double upper = times[i] - source_program.leaf_onset_lag;
      if (!(upper > 0.0)) {
        continue;
      }
      if (need_cdf) {
        scratch.cdf_totals[i] = 0.0;
      }
      const auto lane = static_cast<std::size_t>(
          compiled_frame_lane(lanes, i));
      const auto source_lane = workspace->source_lane(frame, lane);
      const auto &onset_program =
          program.source_programs[
              static_cast<std::size_t>(
                  source_program.onset_source_program_id)];
      const double exact = workspace->source_state().sequence_exact_time(
          source_lane, onset_program.source_id);
      if (std::isfinite(exact)) {
        scratch.lanes[child_count] = compiled_frame_lane(lanes, i);
        scratch.positions[child_count] = i;
        scratch.times[child_count] = exact;
        scratch.parent_times[child_count] = times[i];
        scratch.weights[child_count] = 1.0;
        scratch.exact_flags[child_count] = 1U;
        ++child_count;
        continue;
      }
      const auto nodes = quadrature::map_rule_to_finite_interval<
          quadrature::kDefaultFiniteOrder>(0.0, upper);
      for (std::size_t q = 0; q < quadrature::kDefaultFiniteOrder; ++q) {
        scratch.lanes[child_count] = compiled_frame_lane(lanes, i);
        scratch.positions[child_count] = i;
        scratch.times[child_count] = nodes.nodes[q];
        scratch.parent_times[child_count] = times[i];
        scratch.weights[child_count] = nodes.weights[q];
        scratch.exact_flags[child_count] = 0U;
        ++child_count;
      }
    }
    if (child_count == 0U) {
      continue;
    }

    const auto &onset = evaluate_source_program_lanes(
        program,
        source_program.onset_source_program_id,
        workspace,
        frame,
        scratch.lanes.data(),
        scratch.times.data(),
        child_count,
        kLeafChannelPdf,
        source_workspace,
        depth + 1U);
    scratch.shifted.resize(child_count, shifted_mask);
    source_lane_leaf_program(
        source_program,
        workspace,
        frame,
        scratch.lanes.data(),
        child_count,
        [&](const std::size_t child,
            const ExactLaneLeafBatchView &leaf,
            const int row) {
          return scratch.parent_times[child] - scratch.times[child] -
                 source_program.leaf_onset_lag -
                 leaf.t0(row);
        },
        &scratch.leaf_batch,
        &scratch.shifted);
    for (std::size_t child = 0; child < child_count; ++child) {
      const double onset_density = scratch.exact_flags[child] != 0U
                                       ? 1.0
                                       : onset.pdf[child];
      if (!(onset_density > 0.0)) {
        continue;
      }
      const auto parent = scratch.positions[child];
      const double weight = scratch.weights[child] * onset_density;
      if (need_pdf) {
        scratch.out.pdf[parent] += weight * scratch.shifted.pdf[child];
      }
      if (need_cdf) {
        scratch.cdf_totals[parent] += weight * scratch.shifted.cdf[child];
      }
    }
  }
  for (std::size_t i = 0; i < count; ++i) {
    if (need_pdf) {
      scratch.out.pdf[i] = safe_density(scratch.out.pdf[i]);
    }
    if (need_cdf && std::isfinite(scratch.cdf_totals[i])) {
      const double cdf = clamp_probability(scratch.cdf_totals[i]);
      if ((fill_mask & kLeafChannelCdf) != 0U) {
        scratch.out.cdf[i] = cdf;
      }
      if ((fill_mask & kLeafChannelSurvival) != 0U) {
        scratch.out.survival[i] = clamp_probability(1.0 - cdf);
      }
    }
  }
  return scratch.out;
}

inline SourceLaneFill &evaluate_pool_source_program_lanes(
    const CompiledMathProgram &program,
    const CompiledMathSourceProductProgram &source_program,
    CompiledLaneWorkspace *workspace,
    const CompiledLaneFrame &frame,
    const semantic::Index *lanes,
    const double *times,
    const std::size_t count,
    const std::uint8_t fill_mask,
    SourceLaneWorkspace *source_workspace,
    const std::size_t depth) {
  auto &scratch = source_workspace->layer(depth);
  scratch.out.assign(count, fill_mask);
  const auto members = static_cast<std::size_t>(source_program.member_programs.size);
  const auto k = static_cast<std::size_t>(source_program.pool_k);
  const bool need_pdf = (fill_mask & kLeafChannelPdf) != 0U;
  const bool need_cdf =
      (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
  const auto member_mask = static_cast<std::uint8_t>(
      kLeafChannelCdf | kLeafChannelSurvival |
      (need_pdf ? kLeafChannelPdf : 0U));
  scratch.ensure_pool(count, members, need_pdf);

  for (std::size_t member = 0; member < members; ++member) {
    const auto member_id =
        program.source_program_members[
            static_cast<std::size_t>(source_program.member_programs.offset) +
            member];
    const auto &values = evaluate_source_program_lanes(
        program,
        member_id,
        workspace,
        frame,
        lanes,
        times,
        count,
        member_mask,
        source_workspace,
        depth + 1U);
    for (std::size_t lane = 0; lane < count; ++lane) {
      const auto pos = member * count + lane;
      if (need_pdf) {
        scratch.member_pdf[pos] = values.pdf[lane];
      }
      scratch.member_cdf[pos] = values.cdf[lane];
      scratch.member_survival[pos] = values.survival[lane];
    }
  }

  const auto width = members + 1U;
  const auto table_size = width * width;
  std::fill(
      scratch.prefix.begin(), scratch.prefix.begin() + count * table_size, 0.0);
  if (need_pdf) {
    std::fill(
        scratch.suffix.begin(), scratch.suffix.begin() + count * table_size, 0.0);
  }
  const auto cell = [width, count](const std::size_t row,
                                   const std::size_t col,
                                   const std::size_t lane) {
    return (row * width + col) * count + lane;
  };
  for (std::size_t lane = 0; lane < count; ++lane) {
    scratch.prefix[cell(0U, 0U, lane)] = 1.0;
    if (need_pdf) {
      scratch.suffix[cell(members, 0U, lane)] = 1.0;
    }
  }
  for (std::size_t member = 0; member < members; ++member) {
    for (std::size_t successes = 0; successes <= member; ++successes) {
      for (std::size_t lane = 0; lane < count; ++lane) {
        const auto value_pos = member * count + lane;
        const double value =
            scratch.prefix[cell(member, successes, lane)];
        scratch.prefix[cell(member + 1U, successes, lane)] +=
            value * scratch.member_survival[value_pos];
        scratch.prefix[cell(member + 1U, successes + 1U, lane)] +=
            value * scratch.member_cdf[value_pos];
      }
    }
  }
  if (need_pdf) {
    for (std::size_t member = members; member-- > 0U;) {
      const auto remaining = members - member - 1U;
      for (std::size_t successes = 0; successes <= remaining; ++successes) {
        for (std::size_t lane = 0; lane < count; ++lane) {
          const auto value_pos = member * count + lane;
          const double value =
              scratch.suffix[cell(member + 1U, successes, lane)];
          scratch.suffix[cell(member, successes, lane)] +=
              value * scratch.member_survival[value_pos];
          scratch.suffix[cell(member, successes + 1U, lane)] +=
              value * scratch.member_cdf[value_pos];
        }
      }
    }
  }
  if (need_pdf) {
    for (std::size_t member = 0; member < members; ++member) {
      std::fill(scratch.weights.begin(), scratch.weights.end(), 0.0);
      for (std::size_t left = 0; left < k; ++left) {
        const auto right = k - 1U - left;
        if (right > members - member - 1U) {
          continue;
        }
        for (std::size_t lane = 0; lane < count; ++lane) {
          scratch.weights[lane] +=
              scratch.prefix[cell(member, left, lane)] *
              scratch.suffix[cell(member + 1U, right, lane)];
        }
      }
      for (std::size_t lane = 0; lane < count; ++lane) {
        scratch.out.pdf[lane] +=
            scratch.member_pdf[member * count + lane] *
            scratch.weights[lane];
      }
    }
    for (std::size_t lane = 0; lane < count; ++lane) {
      scratch.out.pdf[lane] = safe_density(scratch.out.pdf[lane]);
    }
  }
  if (need_cdf) {
    std::fill(scratch.weights.begin(), scratch.weights.end(), 0.0);
    for (std::size_t successes = 0; successes < k; ++successes) {
      for (std::size_t lane = 0; lane < count; ++lane) {
        scratch.weights[lane] +=
            scratch.prefix[cell(members, successes, lane)];
      }
    }
    for (std::size_t lane = 0; lane < count; ++lane) {
      const double survival = clamp_probability(scratch.weights[lane]);
      if ((fill_mask & kLeafChannelSurvival) != 0U) {
        scratch.out.survival[lane] = survival;
      }
      if ((fill_mask & kLeafChannelCdf) != 0U) {
        scratch.out.cdf[lane] = clamp_probability(1.0 - survival);
      }
    }
  }
  return scratch.out;
}

inline SourceLaneFill &evaluate_source_program_lanes(
    const CompiledMathProgram &program,
    const semantic::Index source_program_id,
    CompiledLaneWorkspace *workspace,
    const CompiledLaneFrame &frame,
    const semantic::Index *lanes,
    const double *times,
    const std::size_t count,
    const std::uint8_t fill_mask,
    SourceLaneWorkspace *source_workspace,
    const std::size_t depth) {
  auto &scratch = source_workspace->layer(depth);
  scratch.out.resize(count, fill_mask);
  if (source_program_id == semantic::kInvalidIndex || count == 0U) {
    scratch.out.assign(count, fill_mask);
    return scratch.out;
  }
  semantic::Index active_program_id = source_program_id;
  if (!frame.has_sequence_history) {
    const auto &initial = program.source_programs[
        static_cast<std::size_t>(source_program_id)];
    active_program_id = (fill_mask & kLeafChannelPdf) != 0U
                            ? initial.initial_with_pdf_program_id
                            : initial.initial_without_pdf_program_id;
    if (active_program_id == semantic::kInvalidIndex) {
      scratch.out.assign(count, fill_mask);
      return scratch.out;
    }
    if (active_program_id == kInitialCertainSourceProgramId) {
      if ((fill_mask & kLeafChannelPdf) != 0U) {
        std::fill(scratch.out.pdf.begin(), scratch.out.pdf.end(), 0.0);
      }
      if ((fill_mask & kLeafChannelCdf) != 0U) {
        std::fill(scratch.out.cdf.begin(), scratch.out.cdf.end(), 1.0);
      }
      if ((fill_mask & kLeafChannelSurvival) != 0U) {
        std::fill(
            scratch.out.survival.begin(), scratch.out.survival.end(), 0.0);
      }
      return scratch.out;
    }
    const auto &active = program.source_programs[
        static_cast<std::size_t>(active_program_id)];
    if (exact_source_program_relation(active) == ExactRelation::At &&
        (fill_mask & kLeafChannelPdf) != 0U) {
      auto &child = evaluate_source_program_lanes(
          program,
          active.child_program_id,
          workspace,
          frame,
          lanes,
          times,
          count,
          kLeafChannelPdf,
          source_workspace,
          depth + 1U);
      if (fill_mask == kLeafChannelPdf) {
        return child;
      }
      std::copy_n(child.pdf.begin(), count, scratch.out.pdf.begin());
      if ((fill_mask & kLeafChannelCdf) != 0U) {
        std::fill(scratch.out.cdf.begin(), scratch.out.cdf.end(), 1.0);
      }
      if ((fill_mask & kLeafChannelSurvival) != 0U) {
        std::fill(
            scratch.out.survival.begin(), scratch.out.survival.end(), 0.0);
      }
      return scratch.out;
    }
  }
  if (frame.has_sequence_history) {
    while (true) {
      const auto &candidate =
          program.source_programs[
              static_cast<std::size_t>(active_program_id)];
      if (candidate.kind !=
              CompiledMathSourceProductProgramKind::ExactGate &&
          candidate.kind !=
              CompiledMathSourceProductProgramKind::Conditioned) {
        break;
      }
      const auto relation =
          exact_source_program_relation(candidate);
      const bool mixed_at =
          relation == ExactRelation::At &&
          (fill_mask & kLeafChannelPdf) != 0U &&
          (fill_mask & (kLeafChannelCdf | kLeafChannelSurvival)) != 0U;
      if (mixed_at) {
        auto &pdf = evaluate_source_program_lanes(
            program,
            active_program_id,
            workspace,
            frame,
            lanes,
            times,
            count,
            kLeafChannelPdf,
            source_workspace,
            depth + 1U);
        std::copy_n(pdf.pdf.begin(), count, scratch.out.pdf.begin());
        if ((fill_mask & kLeafChannelCdf) != 0U) {
          std::fill(scratch.out.cdf.begin(), scratch.out.cdf.end(), 1.0);
        }
        if ((fill_mask & kLeafChannelSurvival) != 0U) {
          std::fill(
              scratch.out.survival.begin(), scratch.out.survival.end(), 0.0);
        }
        return scratch.out;
      }
      if (exact_source_relation_forces_fill(relation, fill_mask)) {
        for (std::size_t i = 0; i < count; ++i) {
          source_lane_fill_forced(&scratch.out, i, relation);
        }
        return scratch.out;
      }
      bool all_sequence_bounds_inactive = true;
      for (std::size_t i = 0; i < count; ++i) {
        const auto lane = static_cast<std::size_t>(
            compiled_frame_lane(lanes, i));
        const auto source_lane = workspace->source_lane(frame, lane);
        if (!workspace->source_state().sequence_bounds_inactive(
                source_lane, candidate.source_id)) {
          all_sequence_bounds_inactive = false;
          break;
        }
      }
      if (!all_sequence_bounds_inactive) {
        break;
      }
      active_program_id = candidate.child_program_id;
    }
  }
  const auto &source_program =
      program.source_programs[
          static_cast<std::size_t>(active_program_id)];
  switch (source_program.kind) {
  case CompiledMathSourceProductProgramKind::ConstantZero:
    scratch.out.assign(count, fill_mask);
    return scratch.out;
  case CompiledMathSourceProductProgramKind::LeafAbsolute:
    source_lane_leaf_program(
        source_program,
        workspace,
        frame,
        lanes,
        count,
        [&](const std::size_t i,
            const ExactLaneLeafBatchView &leaf,
            const int row) {
          return times[i] - leaf.onset(row) - leaf.t0(row);
        },
        &scratch.leaf_batch,
        &scratch.out);
    return scratch.out;
  case CompiledMathSourceProductProgramKind::ExactGate: {
    scratch.out.assign(count, fill_mask);
    scratch.ensure_items(count);
    std::size_t unresolved = 0U;
    for (std::size_t i = 0; i < count; ++i) {
      const auto lane = static_cast<std::size_t>(
          compiled_frame_lane(lanes, i));
      const auto source_lane = workspace->source_lane(frame, lane);
      const double exact = workspace->source_state().sequence_exact_time(
          source_lane, source_program.source_id);
      if (std::isfinite(exact)) {
        if (times[i] >= exact) {
          source_lane_store_fill(
              &scratch.out,
              i,
              exact_source_certain_fill(fill_mask));
        }
      } else {
        scratch.lanes[unresolved] = compiled_frame_lane(lanes, i);
        scratch.positions[unresolved] = i;
        scratch.times[unresolved] = times[i];
        ++unresolved;
      }
    }
    if (unresolved > 0U) {
      const auto &child = evaluate_source_program_lanes(
          program,
          source_program.child_program_id,
          workspace,
          frame,
          scratch.lanes.data(),
          scratch.times.data(),
          unresolved,
          fill_mask,
          source_workspace,
          depth + 1U);
      for (std::size_t j = 0; j < unresolved; ++j) {
        ExactSourceFill fill;
        if ((fill_mask & kLeafChannelPdf) != 0U) {
          fill.pdf = child.pdf[j];
        }
        if ((fill_mask & kLeafChannelCdf) != 0U) {
          fill.cdf = child.cdf[j];
        }
        if ((fill_mask & kLeafChannelSurvival) != 0U) {
          fill.survival = child.survival[j];
        }
        source_lane_store_fill(&scratch.out, scratch.positions[j], fill);
      }
    }
    return scratch.out;
  }
  case CompiledMathSourceProductProgramKind::Conditioned:
    return evaluate_conditioned_source_program_lanes(
        program,
        source_program,
        workspace,
        frame,
        lanes,
        times,
        count,
        fill_mask,
        source_workspace,
        depth);
  case CompiledMathSourceProductProgramKind::OnsetConvolution:
    return evaluate_onset_source_program_lanes(
        program,
        source_program,
        workspace,
        frame,
        lanes,
        times,
        count,
        fill_mask,
        source_workspace,
        depth);
  case CompiledMathSourceProductProgramKind::PoolKOfN:
    return evaluate_pool_source_program_lanes(
        program,
        source_program,
        workspace,
        frame,
        lanes,
        times,
        count,
        fill_mask,
        source_workspace,
        depth);
  }
  return scratch.out;
}

inline const double *source_lane_fill_values(
    const SourceLaneFill &fill,
    const std::uint8_t channel_mask) noexcept {
  if (channel_mask == kLeafChannelPdf) {
    return fill.pdf.data();
  }
  if (channel_mask == kLeafChannelCdf) {
    return fill.cdf.data();
  }
  return channel_mask == kLeafChannelSurvival
             ? fill.survival.data()
             : nullptr;
}

inline std::vector<double> &source_lane_fill_storage(
    SourceLaneFill *fill,
    const std::uint8_t channel_mask) noexcept {
  if (channel_mask == kLeafChannelPdf) {
    return fill->pdf;
  }
  if (channel_mask == kLeafChannelCdf) {
    return fill->cdf;
  }
  return fill->survival;
}

inline SourceLaneFill &evaluate_source_program_cached_lanes(
    const CompiledMathProgram &program,
    const semantic::Index source_program_id,
    CompiledLaneWorkspace *workspace,
    CompiledLaneFrame *frame,
    const semantic::Index *lanes,
    const double *times,
    const std::size_t count,
    const std::uint8_t fill_mask,
    const bool cache_result,
    SourceLaneWorkspace *source_workspace,
    SourceLaneFill *out) {
  if (source_program_id == semantic::kInvalidIndex) {
    out->assign(count, fill_mask);
    return *out;
  }
  if (count == 0U) {
    out->resize(0U, fill_mask);
    return *out;
  }
  const auto cache_slot =
      cache_result
          ? program.source_program_cache_slots[
                static_cast<std::size_t>(source_program_id)]
          : semantic::kInvalidIndex;
  if (cache_slot == semantic::kInvalidIndex) {
    return evaluate_source_program_lanes(
        program,
        source_program_id,
        workspace,
        *frame,
        lanes,
        times,
        count,
        fill_mask,
        source_workspace,
        1U);
  }
  out->resize(count, fill_mask);
  auto &scratch = source_workspace->layer(0U);
  scratch.ensure_items(count);
  std::size_t missing = 0U;
  for (std::size_t i = 0; i < count; ++i) {
    const auto lane = static_cast<std::size_t>(
        compiled_frame_lane(lanes, i));
    const auto pos = frame->source_program_pos(cache_slot, lane);
    const bool cached =
        frame->source_program_epoch[pos] ==
            frame->source_program_current_epoch &&
        (frame->source_program_valid_mask[pos] & fill_mask) == fill_mask;
    if (cached) {
      if ((fill_mask & kLeafChannelPdf) != 0U) {
        out->pdf[i] = frame->source_program_pdf[pos];
      }
      if ((fill_mask & kLeafChannelCdf) != 0U) {
        out->cdf[i] = frame->source_program_cdf[pos];
      }
      if ((fill_mask & kLeafChannelSurvival) != 0U) {
        out->survival[i] = frame->source_program_survival[pos];
      }
      continue;
    }
    scratch.lanes[missing] = compiled_frame_lane(lanes, i);
    scratch.positions[missing] = i;
    scratch.times[missing] = times[i];
    ++missing;
  }
  if (missing == 0U) {
    return *out;
  }
  const auto &fill = evaluate_source_program_lanes(
      program,
      source_program_id,
      workspace,
      *frame,
      scratch.lanes.data(),
      scratch.times.data(),
      missing,
      fill_mask,
      source_workspace,
      1U);
  for (std::size_t j = 0; j < missing; ++j) {
    const auto output_pos = scratch.positions[j];
    const auto lane = static_cast<std::size_t>(scratch.lanes[j]);
    const auto cache_pos = frame->source_program_pos(cache_slot, lane);
    if (frame->source_program_epoch[cache_pos] !=
        frame->source_program_current_epoch) {
      frame->source_program_epoch[cache_pos] =
          frame->source_program_current_epoch;
      frame->source_program_valid_mask[cache_pos] = 0U;
    }
    frame->source_program_valid_mask[cache_pos] |= fill_mask;
    if ((fill_mask & kLeafChannelPdf) != 0U) {
      const double value = fill.pdf[j];
      frame->source_program_pdf[cache_pos] = value;
      out->pdf[output_pos] = value;
    }
    if ((fill_mask & kLeafChannelCdf) != 0U) {
      const double value = fill.cdf[j];
      frame->source_program_cdf[cache_pos] = value;
      out->cdf[output_pos] = value;
    }
    if ((fill_mask & kLeafChannelSurvival) != 0U) {
      const double value = fill.survival[j];
      frame->source_program_survival[cache_pos] = value;
      out->survival[output_pos] = value;
    }
  }
  return *out;
}

} // namespace detail
} // namespace accumulatr::eval
