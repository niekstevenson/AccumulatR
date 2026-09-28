#pragma once

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <vector>

#include "compile/prep_to_semantic.hpp"
#include "compile/simulation_program.hpp"
#include "leaf/random.hpp"

namespace accumulatr::simulation {

using semantic::Index;
using semantic::ExprKind;
constexpr double infinity = std::numeric_limits<double>::infinity();

struct Readout {
  std::vector<std::string> labels;
  std::vector<double> weights;
  bool missing_rt{false};
};

struct Component {
  std::vector<Index> sources;
  std::vector<Index> expressions;
  std::vector<Index> outcomes;
  int readouts{1};
  int weight_param_index{-1};
};

struct Program {
  explicit Program(const Rcpp::List &prep) : model(compile::compile_prep(prep)) {
    compile::normalize_simulation_expressions(model);
    const Rcpp::List observation(prep["observation"]);
    readouts = Rcpp::as<int>(observation["n_outcomes"]);
    const int global_readouts = Rcpp::as<int>(observation["global_n_outcomes"]);
    const Rcpp::List overrides(observation["component_n_outcomes"]);
    const Rcpp::List outcomes(prep["outcomes"]);
    for (std::size_t i = 0; i < model.outcomes.size(); ++i) {
      const Rcpp::List outcome(outcomes[i]);
      const Rcpp::List options(outcome["options"]);
      Readout policy;
      if (options.containsElementNamed("guess") && !Rf_isNull(options["guess"])) {
        const Rcpp::List guess(options["guess"]);
        policy.labels = Rcpp::as<std::vector<std::string>>(guess["labels"]);
        policy.weights = Rcpp::as<std::vector<double>>(guess["weights"]);
        policy.missing_rt = guess.containsElementNamed("rt_policy") &&
            !Rf_isNull(guess["rt_policy"]) &&
            Rcpp::as<std::string>(guess["rt_policy"]) == "na";
      }
      readout.push_back(std::move(policy));
    }
    for (const auto &definition : model.components) {
      Component component;
      if (!definition.weight_name.empty()) component.weight_param_index = weight_param_count++;
      component.readouts = overrides.containsElementNamed(definition.id.c_str())
          ? Rcpp::as<int>(overrides[definition.id]) : global_readouts;
      std::vector<bool> active(model.leaves.size(), false);
      for (const auto i : definition.active_leaf_indices) active[i] = true;
      std::vector<bool> seen_source(source_count(), false);
      std::vector<bool> seen_expr(model.expr_nodes.size(), false);
      std::function<void(Index)> source = [&](const Index i) {
        if (seen_source[i]) return;
        seen_source[i] = true;
        if (i < static_cast<Index>(model.leaves.size())) {
          if (!active[i]) return;
          const auto &onset = model.leaves[i].onset;
          if (onset.kind != semantic::OnsetKind::Absolute) source(source_id(onset.source));
        } else {
          for (const auto &member : model.pools[i - model.leaves.size()].members)
            source(source_id(member));
        }
        component.sources.push_back(i);
      };
      std::function<void(Index)> expression = [&](const Index i) {
        if (seen_expr[i]) return;
        seen_expr[i] = true;
        const auto &node = model.expr_nodes[i];
        if (node.kind == ExprKind::Event) source(source_id(node.source));
        for (const auto child : node.children) expression(child);
        if (node.kind == ExprKind::Guard) {
          expression(node.reference_child);
          expression(node.blocker_child);
        }
        component.expressions.push_back(i);
      };
      for (std::size_t i = 0; i < model.outcomes.size(); ++i) {
        const auto &outcome = model.outcomes[i];
        const auto &ids = outcome.component_ids;
        if (!ids.empty() && std::find(ids.begin(), ids.end(), definition.id) == ids.end()) continue;
        expression(outcome.expr_root);
        component.outcomes.push_back(i);
      }
      components.push_back(std::move(component));
    }
    for (const auto &trigger : model.triggers)
      trigger_rows.push_back(*std::min_element(trigger.leaf_indices.begin(), trigger.leaf_indices.end()));
    for (const auto &pool : model.pools) max_pool_members = std::max(max_pool_members, pool.members.size());
  }

  Index source_id(const semantic::SourceRef source) const {
    return source.kind == semantic::SourceKind::Leaf
        ? source.index : static_cast<Index>(model.leaves.size()) + source.index;
  }
  std::size_t source_count() const { return model.leaves.size() + model.pools.size(); }

  semantic::SemanticModel model;
  std::vector<Component> components;
  std::vector<Readout> readout;
  std::vector<Index> trigger_rows;
  std::size_t max_pool_members{0};
  int readouts{1};
  int weight_param_count{0};
};

// Readiness is the completion of the prerequisites before the releasing event.
// It resolves positive-mass ties at shared gates without storing source sets.
struct Completion {
  double time{infinity};
  double ready{-infinity};
};

inline bool earlier(const Completion a, const Completion b) {
  return a.time < b.time || (a.time == b.time && a.ready < b.ready);
}

inline int choose(const std::vector<double> &weights) {
  const double total = std::accumulate(weights.begin(), weights.end(), 0.0);
  double draw = R::runif(0.0, total);
  for (std::size_t i = 0; i < weights.size(); ++i) {
    draw -= weights[i];
    if (draw < 0.0) return i;
  }
  return weights.size() - 1;
}

struct Trial {
  explicit Trial(const Program &program)
      : program(program), sources(program.source_count()),
        expressions(program.model.expr_nodes.size()), triggers(program.model.triggers.size()),
        pool_times(program.max_pool_members), weights(program.components.size()) {
    candidates.reserve(program.model.outcomes.size());
  }

  void run(const Component &component, const double *parameters, const int stride,
           const double *onsets) {
    std::fill(sources.begin(), sources.end(), Completion{});
    std::fill(triggers.begin(), triggers.end(), -1);
    const auto &model = program.model;
    for (const auto i : component.sources) {
      auto &value = sources[i];
      if (i < static_cast<Index>(model.leaves.size())) {
        const auto &leaf = model.leaves[i];
        const auto &onset = leaf.onset;
        const double offset = onsets == nullptr ? NA_REAL : onsets[i];
        double start = onset.kind == semantic::OnsetKind::Absolute
            ? (std::isnan(offset) ? onset.absolute_value : offset)
            : sources[program.source_id(onset.source)].time + onset.lag +
                (std::isnan(offset) ? 0.0 : offset);
        if (!std::isfinite(start)) continue;
        const auto trigger = leaf.trigger_index;
        const auto q_row = trigger < 0 ? i : program.trigger_rows[trigger];
        const double q = parameters[q_row];
        bool failed;
        if (trigger < 0) {
          failed = R::runif(0.0, 1.0) < q;
        } else {
          if (triggers[trigger] < 0) triggers[trigger] = R::runif(0.0, 1.0) < q;
          failed = triggers[trigger];
        }
        if (failed) continue;
        const double t0 = parameters[i + stride];
        value.time = start + t0 + leaf::sample_time(leaf.dist, parameters + i + 2 * stride, stride);
      } else {
        const auto &pool = model.pools[i - model.leaves.size()];
        const auto end = pool_times.begin() + pool.members.size();
        const auto kth = pool_times.begin() + pool.k - 1;
        for (std::size_t j = 0; j < pool.members.size(); ++j)
          pool_times[j] = sources[program.source_id(pool.members[j])].time;
        std::nth_element(pool_times.begin(), kth, end);
        value.time = *kth;
        for (std::size_t j = 0; j < pool.members.size(); ++j) {
          const auto member = sources[program.source_id(pool.members[j])];
          pool_times[j] = member.time < value.time ? member.time
              : member.time == value.time ? member.ready : infinity;
        }
        std::nth_element(pool_times.begin(), kth, end);
        value.ready = *kth;
      }
    }
    for (const auto i : component.expressions) {
      const auto &node = model.expr_nodes[i];
      Completion value;
      switch (node.kind) {
      case ExprKind::Event:
        value = sources[program.source_id(node.source)];
        break;
      case ExprKind::And:
        value.time = -infinity;
        for (const auto child : node.children) value.time = std::max(value.time, expressions[child].time);
        for (const auto child : node.children) {
          const auto x = expressions[child];
          value.ready = std::max(value.ready, x.time < value.time ? x.time : x.ready);
        }
        break;
      case ExprKind::Or:
        for (const auto child : node.children)
          if (earlier(expressions[child], value)) value = expressions[child];
        break;
      case ExprKind::Guard:
        value = expressions[node.reference_child];
        if (expressions[node.blocker_child].time < value.time) value.time = infinity;
        break;
      case ExprKind::Not:
        value.time = std::isfinite(expressions[node.children.front()].time) ? infinity : 0.0;
        break;
      case ExprKind::TrueExpr:
        value.time = 0.0;
        break;
      case ExprKind::Impossible:
        break;
      }
      expressions[i] = value;
    }
    candidates.clear();
    for (const auto i : component.outcomes)
      if (std::isfinite(result(i).time)) candidates.push_back(i);
  }

  Completion result(const Index outcome) const {
    return expressions[program.model.outcomes[outcome].expr_root];
  }

  // Full source ancestry is needed only for requested diagnostic output.
  void core(const Index source, Rcpp::NumericVector &times) const {
    const auto &model = program.model;
    if (source < static_cast<Index>(model.leaves.size())) {
      times.push_back(sources[source].time, model.leaves[source].id);
    } else {
      for (const auto member : model.pools[source - model.leaves.size()].members)
        core(program.source_id(member), times);
    }
  }

  Rcpp::List detail(const int component) const {
    const auto &model = program.model;
    Rcpp::List acc_times, pool_times, event_times;
    for (std::size_t i = 0; i < sources.size(); ++i) {
      const bool is_leaf = i < model.leaves.size();
      const auto &id = is_leaf ? model.leaves[i].id : model.pools[i - model.leaves.size()].id;
      Rcpp::NumericVector leaf_times;
      core(i, leaf_times);
      auto event = Rcpp::List::create(Rcpp::Named("time") = sources[i].time,
                                      Rcpp::Named("core") = leaf_times);
      event_times.push_back(event, id);
      if (is_leaf) acc_times.push_back(sources[i].time, id);
      else pool_times.push_back(event, id);
    }
    Rcpp::CharacterVector labels(candidates.size());
    Rcpp::NumericVector times(candidates.size());
    for (std::size_t i = 0; i < candidates.size(); ++i) {
      labels[i] = model.outcomes[candidates[i]].label;
      times[i] = result(candidates[i]).time;
    }
    Rcpp::List outcome_candidates = Rcpp::List::create(
        Rcpp::Named("label") = labels, Rcpp::Named("time") = times);
    outcome_candidates.attr("class") = "data.frame";
    outcome_candidates.attr("row.names") = Rcpp::IntegerVector::create(
        NA_INTEGER, -static_cast<int>(candidates.size()));
    return Rcpp::List::create(
        Rcpp::Named("component") = model.components[component].id,
        Rcpp::Named("acc_times") = acc_times,
        Rcpp::Named("pool_times") = pool_times,
        Rcpp::Named("event_times") = event_times,
        Rcpp::Named("outcome_candidates") = outcome_candidates);
  }

  const Program &program;
  std::vector<Completion> sources, expressions;
  std::vector<int> triggers;
  std::vector<double> pool_times, weights;
  std::vector<Index> candidates;
};

} // namespace accumulatr::simulation
