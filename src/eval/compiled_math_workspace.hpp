#pragma once

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <vector>

#include "compiled_math_types.hpp"

namespace accumulatr::eval {
namespace detail {

struct CompiledMathWorkspace {
  CompiledMathWorkspace() = default;

  explicit CompiledMathWorkspace(const CompiledMathProgram &program) {
    resize(program);
  }

  void resize(const CompiledMathProgram &program) {
    values.assign(program.nodes.size(), 0.0);
    cache_valid.assign(program.nodes.size(), 0U);
    cache_condition_ids.assign(program.nodes.size(), 0);
    cache_times.assign(program.nodes.size(), 0.0);
    cache_evaluators.assign(program.nodes.size(), nullptr);
    cache_values.assign(program.nodes.size(), 0.0);
    const auto source_product_program_count =
        program.integral_kernel_source_product_programs.size();
    source_product_program_epoch.assign(source_product_program_count, 0U);
    source_product_program_valid_mask.assign(
        source_product_program_count, 0U);
    source_product_program_pdf.assign(source_product_program_count, 0.0);
    source_product_program_cdf.assign(source_product_program_count, 0.0);
    source_product_program_survival.assign(source_product_program_count, 1.0);
    source_product_scratch.assign(
        static_cast<std::size_t>(
            program.integral_kernel_source_product_scratch_size),
        0.0);
    time_values.assign(4U, 0.0);
    time_valid.assign(4U, 0U);
    time_valid[static_cast<std::size_t>(CompiledMathTimeSlot::Zero)] = 1U;
  }

  void ensure_size(const CompiledMathProgram &program) {
    if (values.size() < program.nodes.size()) {
      const auto size = program.nodes.size();
      values.resize(size, 0.0);
      cache_valid.resize(size, 0U);
      cache_condition_ids.resize(size, 0);
      cache_times.resize(size, 0.0);
      cache_evaluators.resize(size, nullptr);
      cache_values.resize(size, 0.0);
    }
    if (source_product_program_epoch.size() <
        program.integral_kernel_source_product_programs.size()) {
      const auto size = program.integral_kernel_source_product_programs.size();
      source_product_program_epoch.resize(size, 0U);
      source_product_program_valid_mask.resize(size, 0U);
      source_product_program_pdf.resize(size, 0.0);
      source_product_program_cdf.resize(size, 0.0);
      source_product_program_survival.resize(size, 1.0);
    }
    if (source_product_scratch.size() <
        static_cast<std::size_t>(
            program.integral_kernel_source_product_scratch_size)) {
      source_product_scratch.resize(
          static_cast<std::size_t>(
              program.integral_kernel_source_product_scratch_size),
          0.0);
    }
  }

  void reset_cache() {
    std::fill(cache_valid.begin(), cache_valid.end(), 0U);
    reset_source_product_program_cache();
  }

  void reset_source_product_program_cache() {
    ++source_product_program_current_epoch;
    if (source_product_program_current_epoch == 0U) {
      source_product_program_current_epoch = 1U;
      std::fill(
          source_product_program_epoch.begin(),
          source_product_program_epoch.end(),
          0U);
      std::fill(
          source_product_program_valid_mask.begin(),
          source_product_program_valid_mask.end(),
          0U);
    }
  }

  void set_time(const semantic::Index time_id, const double value) {
    const auto pos = static_cast<std::size_t>(time_id);
    if (time_values.size() <= pos) {
      time_values.resize(pos + 1U, 0.0);
      time_valid.resize(pos + 1U, 0U);
    }
    time_values[pos] = value;
    time_valid[pos] = 1U;
  }

  bool has_time(const semantic::Index time_id) const {
    const auto pos = static_cast<std::size_t>(time_id);
    return pos < time_valid.size() && time_valid[pos] != 0U;
  }

  double time(const semantic::Index time_id) const {
    const auto pos = static_cast<std::size_t>(time_id);
    if (pos >= time_valid.size() || time_valid[pos] == 0U) {
      throw std::runtime_error(
          "compiled math time slot is unbound: " +
          std::to_string(time_id));
    }
    return time_values[pos];
  }

  void copy_times_from(const CompiledMathWorkspace &other) {
    time_values = other.time_values;
    time_valid = other.time_valid;
    used_outcomes = other.used_outcomes;
  }

  class TimeBinding {
  public:
    TimeBinding(CompiledMathWorkspace *workspace,
                const semantic::Index time_id,
                const double value)
        : workspace_(workspace), time_id_(time_id) {
      if (workspace_ == nullptr) {
        return;
      }
      const auto pos = static_cast<std::size_t>(time_id_);
      if (workspace_->time_values.size() <= pos) {
        workspace_->time_values.resize(pos + 1U, 0.0);
        workspace_->time_valid.resize(pos + 1U, 0U);
      }
      had_previous_ = workspace_->time_valid[pos] != 0U;
      previous_ = workspace_->time_values[pos];
      workspace_->time_values[pos] = value;
      workspace_->time_valid[pos] = 1U;
    }

    TimeBinding(const TimeBinding &) = delete;
    TimeBinding &operator=(const TimeBinding &) = delete;

    ~TimeBinding() {
      if (workspace_ == nullptr) {
        return;
      }
      const auto pos = static_cast<std::size_t>(time_id_);
      if (had_previous_) {
        workspace_->time_values[pos] = previous_;
        workspace_->time_valid[pos] = 1U;
      } else {
        workspace_->time_values[pos] = 0.0;
        workspace_->time_valid[pos] = 0U;
      }
    }

  private:
    CompiledMathWorkspace *workspace_{nullptr};
    semantic::Index time_id_{0};
    double previous_{0.0};
    bool had_previous_{false};
  };

  class RebindableTimeBinding {
  public:
    RebindableTimeBinding(CompiledMathWorkspace *workspace,
                          const semantic::Index time_id)
        : workspace_(workspace), time_id_(time_id) {
      if (workspace_ == nullptr) {
        return;
      }
      pos_ = static_cast<std::size_t>(time_id_);
      if (workspace_->time_values.size() <= pos_) {
        workspace_->time_values.resize(pos_ + 1U, 0.0);
        workspace_->time_valid.resize(pos_ + 1U, 0U);
      }
      had_previous_ = workspace_->time_valid[pos_] != 0U;
      previous_ = workspace_->time_values[pos_];
    }

    RebindableTimeBinding(const RebindableTimeBinding &) = delete;
    RebindableTimeBinding &operator=(const RebindableTimeBinding &) = delete;

    ~RebindableTimeBinding() {
      if (workspace_ == nullptr) {
        return;
      }
      if (had_previous_) {
        workspace_->time_values[pos_] = previous_;
        workspace_->time_valid[pos_] = 1U;
      } else {
        workspace_->time_values[pos_] = 0.0;
        workspace_->time_valid[pos_] = 0U;
      }
    }

    void set(const double value) {
      if (workspace_ == nullptr) {
        return;
      }
      workspace_->time_values[pos_] = value;
      workspace_->time_valid[pos_] = 1U;
    }

  private:
    CompiledMathWorkspace *workspace_{nullptr};
    semantic::Index time_id_{0};
    std::size_t pos_{0};
    double previous_{0.0};
    bool had_previous_{false};
  };

  CompiledMathWorkspace &integral_workspace_for(
      const CompiledMathProgram &program) {
    if (!integral_workspace) {
      integral_workspace = std::make_unique<CompiledMathWorkspace>(program);
    } else {
      integral_workspace->reset_cache();
    }
    integral_workspace->copy_times_from(*this);
    return *integral_workspace;
  }

  std::vector<double> values;
  std::vector<std::uint8_t> cache_valid;
  std::vector<std::size_t> cache_condition_ids;
  std::vector<double> cache_times;
  std::vector<const void *> cache_evaluators;
  std::vector<double> cache_values;
  std::vector<std::uint32_t> source_product_program_epoch;
  std::vector<std::uint8_t> source_product_program_valid_mask;
  std::vector<double> source_product_program_pdf;
  std::vector<double> source_product_program_cdf;
  std::vector<double> source_product_program_survival;
  std::uint32_t source_product_program_current_epoch{1U};
  std::vector<double> source_product_scratch;
  std::vector<double> time_values;
  std::vector<std::uint8_t> time_valid;
  std::vector<std::uint8_t> integral_term_open;
  const std::vector<std::uint8_t> *used_outcomes{nullptr};
  std::unique_ptr<CompiledMathWorkspace> integral_workspace;

};

} // namespace detail
} // namespace accumulatr::eval
