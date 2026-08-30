#pragma once

#include <Rcpp.h>

#include <cmath>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>

#include "../semantic/model.hpp"

namespace accumulatr::compile {
namespace detail {

inline std::string as_string(SEXP x) {
  return Rcpp::as<std::string>(x);
}

inline double as_double(SEXP x) {
  return Rcpp::as<double>(x);
}

inline int as_int(SEXP x) {
  return Rcpp::as<int>(x);
}

inline std::vector<std::string> as_string_vector(SEXP x) {
  if (Rf_isNull(x)) {
    return {};
  }
  Rcpp::CharacterVector chr(x);
  std::vector<std::string> out;
  out.reserve(chr.size());
  for (R_xlen_t i = 0; i < chr.size(); ++i) {
    if (chr[i] == NA_STRING) {
      throw std::runtime_error("model ids must not be missing");
    }
    out.push_back(Rcpp::as<std::string>(chr[i]));
  }
  return out;
}

inline int expr_likelihood_id(const Rcpp::RObject &expr_obj) {
  SEXP id_attr = Rf_getAttrib(expr_obj, Rf_install(".lik_id"));
  if (Rf_isNull(id_attr) || Rf_length(id_attr) == 0) {
    return 0;
  }
  if (TYPEOF(id_attr) == INTSXP) {
    const int value = INTEGER(id_attr)[0];
    return value == NA_INTEGER ? 0 : value;
  }
  if (TYPEOF(id_attr) == REALSXP) {
    const double value = REAL(id_attr)[0];
    if (!std::isfinite(value) || value <= 0.0) {
      return 0;
    }
    return static_cast<int>(value);
  }
  return 0;
}

inline std::vector<std::string> internal_param_keys(const std::string &leaf_id,
                                                    leaf::DistKind dist) {
  std::vector<std::string> suffixes;
  switch (dist) {
  case leaf::DistKind::Lognormal:
    suffixes = {"m", "s"};
    break;
  case leaf::DistKind::Gamma:
    suffixes = {"shape", "rate"};
    break;
  case leaf::DistKind::Exgauss:
    suffixes = {"mu", "sigma", "tau"};
    break;
  case leaf::DistKind::LBA:
    suffixes = {"v", "B", "A", "sv"};
    break;
  case leaf::DistKind::RDM:
    suffixes = {"v", "B", "A", "s"};
    break;
  }

  for (auto &suffix : suffixes) {
    suffix = leaf_id + "." + suffix;
  }
  return suffixes;
}

inline semantic::SourceRef source_ref_from_name(
    std::string_view name,
    const std::unordered_map<std::string, semantic::Index> &leaf_index,
    const std::unordered_map<std::string, semantic::Index> &pool_index) {
  auto leaf_it = leaf_index.find(std::string(name));
  if (leaf_it != leaf_index.end()) {
    return semantic::SourceRef{
        semantic::SourceKind::Leaf, leaf_it->second, std::string()};
  }
  auto pool_it = pool_index.find(std::string(name));
  if (pool_it != pool_index.end()) {
    return semantic::SourceRef{
        semantic::SourceKind::Pool, pool_it->second, std::string()};
  }
  throw std::runtime_error("unknown event source '" + std::string(name) + "'");
}

inline semantic::OnsetSpec compile_onset(
    const Rcpp::List &acc,
    const std::unordered_map<std::string, semantic::Index> &leaf_index,
    const std::unordered_map<std::string, semantic::Index> &pool_index) {
  semantic::OnsetSpec out;
  Rcpp::List spec(acc["onset_spec"]);
  const std::string kind = as_string(spec["kind"]);
  if (kind == "absolute") {
    out.kind = semantic::OnsetKind::Absolute;
    out.absolute_value = as_double(spec["value"]);
    return out;
  }
  if (kind != "after") {
    throw std::runtime_error("unsupported onset kind '" + kind + "'");
  }

  const std::string source_name = as_string(spec["source"]);
  const std::string source_kind = as_string(spec["source_kind"]);
  out.lag = as_double(spec["lag"]);
  if (source_kind == "accumulator") {
    out.kind = semantic::OnsetKind::AfterLeaf;
  } else if (source_kind == "pool") {
    out.kind = semantic::OnsetKind::AfterPool;
  } else {
    throw std::runtime_error("unsupported onset source_kind '" + source_kind + "'");
  }
  out.source = source_ref_from_name(source_name, leaf_index, pool_index);
  return out;
}

inline semantic::ExprKind expr_kind_from_string(std::string_view kind) {
  if (kind == "event") return semantic::ExprKind::Event;
  if (kind == "and") return semantic::ExprKind::And;
  if (kind == "or") return semantic::ExprKind::Or;
  if (kind == "not") return semantic::ExprKind::Not;
  if (kind == "guard") return semantic::ExprKind::Guard;
  if (kind == "impossible") return semantic::ExprKind::Impossible;
  if (kind == "true") return semantic::ExprKind::TrueExpr;
  throw std::runtime_error("unsupported expr kind '" + std::string(kind) + "'");
}

inline semantic::Index compile_expr(
    const Rcpp::RObject &expr_obj, semantic::SemanticModel *model,
    const std::unordered_map<std::string, semantic::Index> &leaf_index,
    const std::unordered_map<std::string, semantic::Index> &pool_index,
    std::unordered_map<int, semantic::Index> *expr_id_index = nullptr) {
  if (expr_obj.isNULL()) {
    throw std::runtime_error("expression node must not be NULL");
  }

  const int likelihood_id = expr_likelihood_id(expr_obj);
  if (likelihood_id > 0 && expr_id_index != nullptr) {
    const auto found = expr_id_index->find(likelihood_id);
    if (found != expr_id_index->end()) {
      return found->second;
    }
  }

  Rcpp::List expr(expr_obj);
  semantic::ExprNode node;
  node.kind = expr_kind_from_string(as_string(expr["kind"]));

  switch (node.kind) {
  case semantic::ExprKind::Event: {
    const std::string source_name = as_string(expr["source"]);
    node.source = source_ref_from_name(source_name, leaf_index, pool_index);
    break;
  }
  case semantic::ExprKind::And:
  case semantic::ExprKind::Or: {
    Rcpp::List args(expr["args"]);
    node.children.reserve(args.size());
    for (R_xlen_t i = 0; i < args.size(); ++i) {
      node.children.push_back(
          compile_expr(args[i], model, leaf_index, pool_index, expr_id_index));
    }
    break;
  }
  case semantic::ExprKind::Not: {
    node.children.push_back(
        compile_expr(expr["arg"], model, leaf_index, pool_index, expr_id_index));
    break;
  }
  case semantic::ExprKind::Guard: {
    node.reference_child =
        compile_expr(expr["reference"], model, leaf_index, pool_index, expr_id_index);
    node.blocker_child =
        compile_expr(expr["blocker"], model, leaf_index, pool_index, expr_id_index);
    if (expr.containsElementNamed("unless") && !Rf_isNull(expr["unless"])) {
      Rcpp::List unless_list(expr["unless"]);
      node.unless_children.reserve(unless_list.size());
      for (R_xlen_t i = 0; i < unless_list.size(); ++i) {
        node.unless_children.push_back(
            compile_expr(unless_list[i], model, leaf_index, pool_index, expr_id_index));
      }
    }
    break;
  }
  case semantic::ExprKind::Impossible:
  case semantic::ExprKind::TrueExpr:
    break;
  }

  model->expr_nodes.push_back(std::move(node));
  const auto index = static_cast<semantic::Index>(model->expr_nodes.size() - 1);
  if (likelihood_id > 0 && expr_id_index != nullptr) {
    (*expr_id_index)[likelihood_id] = index;
  }
  return index;
}

} // namespace detail

inline semantic::SemanticModel compile_prep(const Rcpp::List &prep) {
  semantic::SemanticModel model;

  if (!prep.containsElementNamed("accumulators") ||
      !prep.containsElementNamed("outcomes")) {
    throw std::runtime_error("prep must contain accumulators and outcomes");
  }

  Rcpp::List accs(prep["accumulators"]);
  Rcpp::CharacterVector acc_names(accs.names());
  model.leaves.reserve(accs.size());
  for (R_xlen_t i = 0; i < accs.size(); ++i) {
    Rcpp::List acc(accs[i]);
    const std::string leaf_id = Rcpp::as<std::string>(acc_names[i]);
    leaf::DistKind dist{};
    const std::string dist_name = detail::as_string(acc["dist"]);
    if (!leaf::try_parse_dist_kind(dist_name, &dist)) {
      throw std::runtime_error("unknown distribution '" + dist_name + "'");
    }

    semantic::LeafSpec leaf_spec;
    leaf_spec.id = leaf_id;
    leaf_spec.dist = dist;
    leaf_spec.params.dist_param_names = detail::internal_param_keys(leaf_id, dist);
    leaf_spec.params.t0_name = leaf_id + ".t0";
    model.leaves.push_back(std::move(leaf_spec));
  }

  std::unordered_map<std::string, semantic::Index> leaf_index;
  leaf_index.reserve(model.leaves.size());
  for (semantic::Index i = 0; i < static_cast<semantic::Index>(model.leaves.size());
       ++i) {
    leaf_index.emplace(model.leaves[static_cast<std::size_t>(i)].id, i);
  }

  Rcpp::List pools(prep["pools"]);
  Rcpp::CharacterVector pool_names;
  if (pools.size() > 0) {
    pool_names = pools.names();
  }
  model.pools.reserve(pools.size());
  for (R_xlen_t i = 0; i < pools.size(); ++i) {
    Rcpp::List pool(pools[i]);
    semantic::PoolSpec pool_spec;
    pool_spec.id = Rcpp::as<std::string>(pool_names[i]);
    pool_spec.k = detail::as_int(pool["k"]);
    model.pools.push_back(std::move(pool_spec));
  }

  std::unordered_map<std::string, semantic::Index> pool_index;
  pool_index.reserve(model.pools.size());
  for (semantic::Index i = 0; i < static_cast<semantic::Index>(model.pools.size());
       ++i) {
    pool_index.emplace(model.pools[static_cast<std::size_t>(i)].id, i);
  }

  for (R_xlen_t i = 0; i < accs.size(); ++i) {
    Rcpp::List acc(accs[i]);
    model.leaves[static_cast<std::size_t>(i)].onset =
        detail::compile_onset(acc, leaf_index, pool_index);
  }

  for (R_xlen_t i = 0; i < pools.size(); ++i) {
    Rcpp::List pool(pools[i]);
    const auto member_names = detail::as_string_vector(pool["members"]);
    auto &members = model.pools[static_cast<std::size_t>(i)].members;
    members.reserve(member_names.size());
    for (const auto &name : member_names) {
      members.push_back(detail::source_ref_from_name(name, leaf_index, pool_index));
    }
  }

  Rcpp::List triggers(prep["shared_triggers"]);
  Rcpp::CharacterVector trigger_names;
  if (triggers.size() > 0) {
    trigger_names = triggers.names();
  }
  model.triggers.reserve(triggers.size());
  for (R_xlen_t i = 0; i < triggers.size(); ++i) {
    Rcpp::List trigger(triggers[i]);
    semantic::TriggerSpec trigger_spec;
    trigger_spec.id = Rcpp::as<std::string>(trigger_names[i]);
    const auto member_names = detail::as_string_vector(trigger["members"]);
    for (const auto &name : member_names) {
      auto it = leaf_index.find(name);
      if (it == leaf_index.end()) {
        throw std::runtime_error("shared trigger member '" + name +
                                 "' is not a known leaf");
      }
      trigger_spec.leaf_indices.push_back(it->second);
    }
    model.triggers.push_back(std::move(trigger_spec));
  }

  for (semantic::Index i = 0; i < static_cast<semantic::Index>(model.triggers.size());
       ++i) {
    for (const auto leaf_i :
         model.triggers[static_cast<std::size_t>(i)].leaf_indices) {
      model.leaves[static_cast<std::size_t>(leaf_i)].trigger_index = i;
    }
  }

  Rcpp::List outcomes(prep["outcomes"]);
  model.outcomes.reserve(outcomes.size());
  std::unordered_map<int, semantic::Index> expr_id_index;
  Rcpp::CharacterVector outcome_names(outcomes.names());
  for (R_xlen_t i = 0; i < outcomes.size(); ++i) {
    Rcpp::List outcome(outcomes[i]);
    semantic::OutcomeSpec outcome_spec;
    outcome_spec.label = Rcpp::as<std::string>(outcome_names[i]);
    outcome_spec.expr_root = detail::compile_expr(
        outcome["expr"], &model, leaf_index, pool_index, &expr_id_index);
    Rcpp::List options(outcome["options"]);
    if (options.containsElementNamed("component") &&
        !Rf_isNull(options["component"])) {
      outcome_spec.component_ids =
          detail::as_string_vector(options["component"]);
    }
    if (options.containsElementNamed("map_outcome_to") &&
        !Rf_isNull(options["map_outcome_to"])) {
      Rcpp::CharacterVector target(options["map_outcome_to"]);
      if (target[0] == NA_STRING) {
        outcome_spec.mapping.maps_to_missing = true;
      } else {
        outcome_spec.mapping.observed_label = Rcpp::as<std::string>(target[0]);
      }
    }
    outcome_spec.has_guess =
        options.containsElementNamed("guess") &&
        !Rf_isNull(options["guess"]);
    model.outcomes.push_back(std::move(outcome_spec));
  }

  Rcpp::List components(prep["components"]);
  model.component_mode = detail::as_string(components["mode"]);
  model.component_reference = detail::as_string(components["reference"]);

  const auto component_ids = detail::as_string_vector(components["ids"]);
  Rcpp::NumericVector weights(components["weights"]);
  Rcpp::List attrs(components["attrs"]);

  std::unordered_map<std::string, std::vector<semantic::Index>> active_map;
  for (semantic::Index i = 0;
       i < static_cast<semantic::Index>(model.leaves.size()); ++i) {
    const Rcpp::List acc(accs[static_cast<R_xlen_t>(i)]);
    const auto active_ids = detail::as_string_vector(acc["components"]);
    for (const auto &component_id : active_ids) {
      active_map[component_id].push_back(i);
    }
  }

  model.components.reserve(component_ids.size());
  for (std::size_t i = 0; i < component_ids.size(); ++i) {
    semantic::ComponentSpec component;
    component.id = component_ids[i];
    component.weight = weights[static_cast<R_xlen_t>(i)];
    if (component.id == "__default__") {
      component.active_leaf_indices.resize(model.leaves.size());
      for (std::size_t j = 0; j < model.leaves.size(); ++j) {
        component.active_leaf_indices[j] = static_cast<semantic::Index>(j);
      }
    } else {
      component.active_leaf_indices = active_map[component.id];
    }
    Rcpp::List attr(attrs[component.id]);
    if (attr.containsElementNamed("weight_param")) {
      component.weight_name = detail::as_string(attr["weight_param"]);
    }
    if (attr.containsElementNamed("n_outcomes") &&
        !Rf_isNull(attr["n_outcomes"])) {
      component.n_outcomes_override = detail::as_int(attr["n_outcomes"]);
    }
    model.components.push_back(std::move(component));
  }

  Rcpp::List observation(prep["observation"]);
  model.observation.mode = semantic::ObservationMode::TopK;
  model.observation.n_outcomes = detail::as_int(observation["n_outcomes"]);
  model.observation.global_n_outcomes =
      detail::as_int(observation["global_n_outcomes"]);

  return model;
}

} // namespace accumulatr::compile
