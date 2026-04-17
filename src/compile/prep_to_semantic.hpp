#pragma once

#include <Rcpp.h>

#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>

#include "../semantic/model.hpp"

namespace accumulatr::compile {
namespace detail {

inline std::string as_string(SEXP x, const std::string &default_value = "") {
  if (Rf_isNull(x)) {
    return default_value;
  }
  if (TYPEOF(x) == STRSXP) {
    Rcpp::CharacterVector chr(x);
    if (chr.size() == 0 || chr[0] == NA_STRING) {
      return default_value;
    }
    return Rcpp::as<std::string>(chr[0]);
  }
  try {
    return Rcpp::as<std::string>(x);
  } catch (...) {
    return default_value;
  }
}

inline double as_double(SEXP x, double default_value = 0.0) {
  if (Rf_isNull(x)) {
    return default_value;
  }
  try {
    return Rcpp::as<double>(x);
  } catch (...) {
    return default_value;
  }
}

inline int as_int(SEXP x, int default_value = 0) {
  if (Rf_isNull(x)) {
    return default_value;
  }
  try {
    return Rcpp::as<int>(x);
  } catch (...) {
    return default_value;
  }
}

inline bool as_bool(SEXP x, bool default_value = false) {
  if (Rf_isNull(x)) {
    return default_value;
  }
  try {
    return Rcpp::as<bool>(x);
  } catch (...) {
    return default_value;
  }
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
      continue;
    }
    out.push_back(Rcpp::as<std::string>(chr[i]));
  }
  return out;
}

inline bool has_named_element(const Rcpp::List &list, const char *name) {
  SEXP names_sexp = list.names();
  if (Rf_isNull(names_sexp)) {
    return false;
  }
  Rcpp::CharacterVector names(names_sexp);
  for (R_xlen_t i = 0; i < names.size(); ++i) {
    if (names[i] == NA_STRING) {
      continue;
    }
    if (std::string(Rcpp::as<std::string>(names[i])) == name) {
      return true;
    }
  }
  return false;
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
  return semantic::SourceRef{
      semantic::SourceKind::Special, semantic::kInvalidIndex, std::string(name)};
}

inline semantic::OnsetSpec compile_onset(
    const Rcpp::List &acc,
    const std::unordered_map<std::string, semantic::Index> &leaf_index,
    const std::unordered_map<std::string, semantic::Index> &pool_index) {
  semantic::OnsetSpec out;
  if (!acc.containsElementNamed("onset_spec") || Rf_isNull(acc["onset_spec"])) {
    out.kind = semantic::OnsetKind::Absolute;
    out.absolute_value = acc.containsElementNamed("onset")
                             ? as_double(acc["onset"], 0.0)
                             : 0.0;
    return out;
  }

  Rcpp::List spec(acc["onset_spec"]);
  const std::string kind = as_string(spec["kind"], "absolute");
  if (kind == "absolute") {
    out.kind = semantic::OnsetKind::Absolute;
    out.absolute_value = as_double(spec["value"], 0.0);
    return out;
  }
  if (kind != "after") {
    throw std::runtime_error("unsupported onset kind '" + kind + "'");
  }

  const std::string source_name = as_string(spec["source"]);
  const std::string source_kind = as_string(spec["source_kind"]);
  out.lag = as_double(spec["lag"], 0.0);
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
    const std::unordered_map<std::string, semantic::Index> &pool_index) {
  if (expr_obj.isNULL()) {
    throw std::runtime_error("expression node must not be NULL");
  }

  Rcpp::List expr(expr_obj);
  semantic::ExprNode node;
  node.kind = expr_kind_from_string(as_string(expr["kind"]));

  switch (node.kind) {
  case semantic::ExprKind::Event: {
    const std::string source_name = as_string(expr["source"]);
    node.source = source_ref_from_name(source_name, leaf_index, pool_index);
    if (expr.containsElementNamed("k") && !Rf_isNull(expr["k"])) {
      node.event_k = as_int(expr["k"], 0);
    }
    break;
  }
  case semantic::ExprKind::And:
  case semantic::ExprKind::Or: {
    Rcpp::List args(expr["args"]);
    node.children.reserve(args.size());
    for (R_xlen_t i = 0; i < args.size(); ++i) {
      node.children.push_back(
          compile_expr(args[i], model, leaf_index, pool_index));
    }
    break;
  }
  case semantic::ExprKind::Not: {
    node.children.push_back(
        compile_expr(expr["arg"], model, leaf_index, pool_index));
    break;
  }
  case semantic::ExprKind::Guard: {
    node.reference_child =
        compile_expr(expr["reference"], model, leaf_index, pool_index);
    node.blocker_child =
        compile_expr(expr["blocker"], model, leaf_index, pool_index);
    if (expr.containsElementNamed("unless") && !Rf_isNull(expr["unless"])) {
      Rcpp::List unless_list(expr["unless"]);
      node.unless_children.reserve(unless_list.size());
      for (R_xlen_t i = 0; i < unless_list.size(); ++i) {
        node.unless_children.push_back(
            compile_expr(unless_list[i], model, leaf_index, pool_index));
      }
    }
    break;
  }
  case semantic::ExprKind::Impossible:
  case semantic::ExprKind::TrueExpr:
    break;
  }

  model->expr_nodes.push_back(std::move(node));
  return static_cast<semantic::Index>(model->expr_nodes.size() - 1);
}

inline Rcpp::List to_r_list(const semantic::SemanticModel &model) {
  Rcpp::List leaves(model.leaves.size());
  for (std::size_t i = 0; i < model.leaves.size(); ++i) {
    const auto &leaf = model.leaves[i];
    std::string onset_kind = "absolute";
    std::string onset_source_kind;
    std::string onset_source_id;
    switch (leaf.onset.kind) {
    case semantic::OnsetKind::Absolute:
      onset_kind = "absolute";
      break;
    case semantic::OnsetKind::AfterLeaf:
      onset_kind = "after";
      onset_source_kind = "leaf";
      onset_source_id = model.leaves[static_cast<std::size_t>(leaf.onset.source.index)].id;
      break;
    case semantic::OnsetKind::AfterPool:
      onset_kind = "after";
      onset_source_kind = "pool";
      onset_source_id = model.pools[static_cast<std::size_t>(leaf.onset.source.index)].id;
      break;
    }
    leaves[i] = Rcpp::List::create(
        Rcpp::Named("id") = leaf.id,
        Rcpp::Named("dist") = std::string(leaf::to_string(leaf.dist)),
        Rcpp::Named("trigger_index") =
            leaf.trigger_index == semantic::kInvalidIndex ? NA_INTEGER
                                                          : leaf.trigger_index + 1,
        Rcpp::Named("param_keys") = leaf.params.dist_param_names,
        Rcpp::Named("q_key") = leaf.params.q_name,
        Rcpp::Named("t0_key") = leaf.params.t0_name,
        Rcpp::Named("onset") = Rcpp::List::create(
            Rcpp::Named("kind") = onset_kind,
            Rcpp::Named("absolute_value") = leaf.onset.absolute_value,
            Rcpp::Named("lag") = leaf.onset.lag,
            Rcpp::Named("source_kind") = onset_source_kind,
            Rcpp::Named("source_id") = onset_source_id));
  }

  Rcpp::List pools(model.pools.size());
  for (std::size_t i = 0; i < model.pools.size(); ++i) {
    const auto &pool = model.pools[i];
    Rcpp::CharacterVector members(pool.members.size());
    for (std::size_t j = 0; j < pool.members.size(); ++j) {
      const auto &member = pool.members[j];
      if (member.kind == semantic::SourceKind::Leaf) {
        members[j] = model.leaves[static_cast<std::size_t>(member.index)].id;
      } else if (member.kind == semantic::SourceKind::Pool) {
        members[j] = model.pools[static_cast<std::size_t>(member.index)].id;
      } else {
        members[j] = member.special_id;
      }
    }
    pools[i] = Rcpp::List::create(
        Rcpp::Named("id") = pool.id,
        Rcpp::Named("k") = pool.k,
        Rcpp::Named("members") = members);
  }

  Rcpp::List triggers(model.triggers.size());
  for (std::size_t i = 0; i < model.triggers.size(); ++i) {
    const auto &trigger = model.triggers[i];
    Rcpp::CharacterVector members(trigger.leaf_indices.size());
    for (std::size_t j = 0; j < trigger.leaf_indices.size(); ++j) {
      members[j] = model.leaves[static_cast<std::size_t>(trigger.leaf_indices[j])].id;
    }
    triggers[i] = Rcpp::List::create(
        Rcpp::Named("id") = trigger.id,
        Rcpp::Named("kind") =
            trigger.kind == semantic::TriggerKind::Shared ? "shared"
                                                          : "independent",
        Rcpp::Named("members") = members,
        Rcpp::Named("q_name") = trigger.q_name,
        Rcpp::Named("fixed_q") = trigger.fixed_q,
        Rcpp::Named("has_fixed_q") = trigger.has_fixed_q);
  }

  Rcpp::List expr_nodes(model.expr_nodes.size());
  for (std::size_t i = 0; i < model.expr_nodes.size(); ++i) {
    const auto &node = model.expr_nodes[i];
    std::string kind;
    switch (node.kind) {
    case semantic::ExprKind::Event:
      kind = "event";
      break;
    case semantic::ExprKind::And:
      kind = "and";
      break;
    case semantic::ExprKind::Or:
      kind = "or";
      break;
    case semantic::ExprKind::Not:
      kind = "not";
      break;
    case semantic::ExprKind::Guard:
      kind = "guard";
      break;
    case semantic::ExprKind::Impossible:
      kind = "impossible";
      break;
    case semantic::ExprKind::TrueExpr:
      kind = "true";
      break;
    }

    std::string source_kind;
    std::string source_id;
    if (node.kind == semantic::ExprKind::Event) {
      if (node.source.kind == semantic::SourceKind::Leaf) {
        source_kind = "leaf";
        source_id = model.leaves[static_cast<std::size_t>(node.source.index)].id;
      } else if (node.source.kind == semantic::SourceKind::Pool) {
        source_kind = "pool";
        source_id = model.pools[static_cast<std::size_t>(node.source.index)].id;
      } else {
        source_kind = "special";
        source_id = node.source.special_id;
      }
    }

    Rcpp::IntegerVector children(node.children.size());
    for (std::size_t j = 0; j < node.children.size(); ++j) {
      children[j] = node.children[j] + 1;
    }
    Rcpp::IntegerVector unless_children(node.unless_children.size());
    for (std::size_t j = 0; j < node.unless_children.size(); ++j) {
      unless_children[j] = node.unless_children[j] + 1;
    }

    expr_nodes[i] = Rcpp::List::create(
        Rcpp::Named("kind") = kind,
        Rcpp::Named("source_kind") = source_kind,
        Rcpp::Named("source_id") = source_id,
        Rcpp::Named("event_k") = node.event_k == 0 ? NA_INTEGER : node.event_k,
        Rcpp::Named("children") = children,
        Rcpp::Named("reference_child") =
            node.reference_child == semantic::kInvalidIndex
                ? NA_INTEGER
                : node.reference_child + 1,
        Rcpp::Named("blocker_child") =
            node.blocker_child == semantic::kInvalidIndex
                ? NA_INTEGER
                : node.blocker_child + 1,
        Rcpp::Named("unless_children") = unless_children);
  }

  Rcpp::List outcomes(model.outcomes.size());
  for (std::size_t i = 0; i < model.outcomes.size(); ++i) {
    const auto &outcome = model.outcomes[i];
    outcomes[i] = Rcpp::List::create(
        Rcpp::Named("label") = outcome.label,
        Rcpp::Named("expr_root") = outcome.expr_root + 1,
        Rcpp::Named("component_ids") = outcome.component_ids,
        Rcpp::Named("has_guess") = outcome.has_guess,
        Rcpp::Named("outcome_class") = outcome.outcome_class,
        Rcpp::Named("mapping") = Rcpp::List::create(
            Rcpp::Named("maps_to_missing") = outcome.mapping.maps_to_missing,
            Rcpp::Named("observed_label") = outcome.mapping.observed_label));
  }

  Rcpp::List components(model.components.size());
  for (std::size_t i = 0; i < model.components.size(); ++i) {
    const auto &component = model.components[i];
    Rcpp::CharacterVector active_leaves(component.active_leaf_indices.size());
    for (std::size_t j = 0; j < component.active_leaf_indices.size(); ++j) {
      active_leaves[j] =
          model.leaves[static_cast<std::size_t>(component.active_leaf_indices[j])].id;
    }
    components[i] = Rcpp::List::create(
        Rcpp::Named("id") = component.id,
        Rcpp::Named("active_leaves") = active_leaves,
        Rcpp::Named("weight") = component.weight,
        Rcpp::Named("weight_name") = component.weight_name,
        Rcpp::Named("n_outcomes_override") =
            component.n_outcomes_override == 0 ? NA_INTEGER
                                               : component.n_outcomes_override);
  }

  const auto issues = semantic::validate_basic(model);
  Rcpp::CharacterVector issue_messages(issues.size());
  for (std::size_t i = 0; i < issues.size(); ++i) {
    issue_messages[i] = issues[i].message;
  }

  return Rcpp::List::create(
      Rcpp::Named("leaves") = leaves,
      Rcpp::Named("pools") = pools,
      Rcpp::Named("triggers") = triggers,
      Rcpp::Named("expr_nodes") = expr_nodes,
      Rcpp::Named("outcomes") = outcomes,
      Rcpp::Named("components") = components,
      Rcpp::Named("observation") = Rcpp::List::create(
          Rcpp::Named("mode") = "top_k",
          Rcpp::Named("n_outcomes") = model.observation.n_outcomes,
          Rcpp::Named("global_n_outcomes") =
              model.observation.global_n_outcomes),
      Rcpp::Named("component_mode") = model.component_mode,
      Rcpp::Named("component_reference") = model.component_reference,
      Rcpp::Named("validation_issues") = issue_messages);
}

} // namespace detail

inline semantic::SemanticModel compile_prep(const Rcpp::List &prep) {
  semantic::SemanticModel model;

  if (!prep.containsElementNamed("accumulators") ||
      !prep.containsElementNamed("outcomes")) {
    throw std::runtime_error("prep must contain accumulators and outcomes");
  }

  Rcpp::List accs(prep["accumulators"]);
  SEXP acc_names_sexp = accs.names();
  Rcpp::CharacterVector acc_names =
      Rf_isNull(acc_names_sexp) ? Rcpp::CharacterVector()
                                : Rcpp::CharacterVector(acc_names_sexp);
  model.leaves.reserve(accs.size());
  for (R_xlen_t i = 0; i < accs.size(); ++i) {
    Rcpp::List acc(accs[i]);
    const std::string leaf_id =
        !acc_names.isNULL() && acc_names[i] != NA_STRING
            ? Rcpp::as<std::string>(acc_names[i])
            : detail::as_string(acc["id"]);
    leaf::DistKind dist{};
    const std::string dist_name = detail::as_string(acc["dist"]);
    if (!leaf::try_parse_dist_kind(dist_name, &dist)) {
      throw std::runtime_error("unknown distribution '" + dist_name + "'");
    }

    semantic::LeafSpec leaf_spec;
    leaf_spec.id = leaf_id;
    leaf_spec.dist = dist;
    leaf_spec.params.dist_param_names = detail::internal_param_keys(leaf_id, dist);
    leaf_spec.params.q_name = leaf_id + ".q";
    leaf_spec.params.t0_name = leaf_id + ".t0";
    model.leaves.push_back(std::move(leaf_spec));
  }

  std::unordered_map<std::string, semantic::Index> leaf_index;
  leaf_index.reserve(model.leaves.size());
  for (semantic::Index i = 0; i < static_cast<semantic::Index>(model.leaves.size());
       ++i) {
    leaf_index.emplace(model.leaves[static_cast<std::size_t>(i)].id, i);
  }

  if (prep.containsElementNamed("pools") && !Rf_isNull(prep["pools"])) {
    Rcpp::List pools(prep["pools"]);
    SEXP pool_names_sexp = pools.names();
    Rcpp::CharacterVector pool_names =
        Rf_isNull(pool_names_sexp) ? Rcpp::CharacterVector()
                                   : Rcpp::CharacterVector(pool_names_sexp);
    model.pools.reserve(pools.size());
    for (R_xlen_t i = 0; i < pools.size(); ++i) {
      Rcpp::List pool(pools[i]);
      semantic::PoolSpec pool_spec;
      pool_spec.id = !pool_names.isNULL() && pool_names[i] != NA_STRING
                         ? Rcpp::as<std::string>(pool_names[i])
                         : detail::as_string(pool["id"]);
      pool_spec.k = detail::as_int(pool["k"], 1);
      model.pools.push_back(std::move(pool_spec));
    }
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

  if (prep.containsElementNamed("pools") && !Rf_isNull(prep["pools"])) {
    Rcpp::List pools(prep["pools"]);
    for (R_xlen_t i = 0; i < pools.size(); ++i) {
      Rcpp::List pool(pools[i]);
      const auto member_names = detail::as_string_vector(pool["members"]);
      auto &members = model.pools[static_cast<std::size_t>(i)].members;
      members.reserve(member_names.size());
      for (const auto &name : member_names) {
        members.push_back(detail::source_ref_from_name(name, leaf_index, pool_index));
      }
    }
  }

  if (prep.containsElementNamed("shared_triggers") &&
      !Rf_isNull(prep["shared_triggers"])) {
    Rcpp::List triggers(prep["shared_triggers"]);
    SEXP trigger_names_sexp = triggers.names();
    Rcpp::CharacterVector trigger_names =
        Rf_isNull(trigger_names_sexp) ? Rcpp::CharacterVector()
                                      : Rcpp::CharacterVector(trigger_names_sexp);
    model.triggers.reserve(triggers.size());
    for (R_xlen_t i = 0; i < triggers.size(); ++i) {
      Rcpp::List trigger(triggers[i]);
      semantic::TriggerSpec trigger_spec;
      trigger_spec.id =
          !trigger_names.isNULL() && trigger_names[i] != NA_STRING
              ? Rcpp::as<std::string>(trigger_names[i])
              : detail::as_string(trigger["id"]);
      const std::string draw = detail::as_string(trigger["draw"], "shared");
      trigger_spec.kind = draw == "independent" ? semantic::TriggerKind::Independent
                                                : semantic::TriggerKind::Shared;
      trigger_spec.q_name = trigger_spec.id;
      if (trigger.containsElementNamed("q") && !Rf_isNull(trigger["q"])) {
        trigger_spec.fixed_q = detail::as_double(trigger["q"], 0.0);
        trigger_spec.has_fixed_q = true;
      }
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
  }

  for (semantic::Index i = 0; i < static_cast<semantic::Index>(model.triggers.size());
       ++i) {
    for (const auto leaf_i :
         model.triggers[static_cast<std::size_t>(i)].leaf_indices) {
      model.leaves[static_cast<std::size_t>(leaf_i)].trigger_index = i;
      model.leaves[static_cast<std::size_t>(leaf_i)].params.q_name =
          model.triggers[static_cast<std::size_t>(i)].q_name;
    }
  }

  Rcpp::List outcomes(prep["outcomes"]);
  model.outcomes.reserve(outcomes.size());
  SEXP outcome_names_sexp = outcomes.names();
  Rcpp::CharacterVector outcome_names =
      Rf_isNull(outcome_names_sexp) ? Rcpp::CharacterVector()
                                    : Rcpp::CharacterVector(outcome_names_sexp);
  for (R_xlen_t i = 0; i < outcomes.size(); ++i) {
    Rcpp::List outcome(outcomes[i]);
    semantic::OutcomeSpec outcome_spec;
    outcome_spec.label =
        !outcome_names.isNULL() && outcome_names[i] != NA_STRING
            ? Rcpp::as<std::string>(outcome_names[i])
            : detail::as_string(outcome["label"]);
    outcome_spec.expr_root = detail::compile_expr(
        outcome["expr"], &model, leaf_index, pool_index);
    if (outcome.containsElementNamed("options") && !Rf_isNull(outcome["options"])) {
      Rcpp::List options(outcome["options"]);
      if (detail::has_named_element(options, "component") &&
          !Rf_isNull(options["component"])) {
        outcome_spec.component_ids =
            detail::as_string_vector(options["component"]);
      }
      if (detail::has_named_element(options, "map_outcome_to") &&
          !Rf_isNull(options["map_outcome_to"])) {
        SEXP map_val = options["map_outcome_to"];
        if (TYPEOF(map_val) == STRSXP) {
          Rcpp::CharacterVector chr(map_val);
          if (chr.size() > 0 && chr[0] == NA_STRING) {
            outcome_spec.mapping.maps_to_missing = true;
          } else if (chr.size() > 0) {
            outcome_spec.mapping.observed_label = Rcpp::as<std::string>(chr[0]);
          }
        }
      }
      outcome_spec.has_guess =
          detail::has_named_element(options, "guess") &&
          !Rf_isNull(options["guess"]);
      if (detail::has_named_element(options, "class") &&
          !Rf_isNull(options["class"])) {
        outcome_spec.outcome_class = detail::as_string(options["class"]);
      }
    }
    model.outcomes.push_back(std::move(outcome_spec));
  }

  if (prep.containsElementNamed("components") && !Rf_isNull(prep["components"])) {
    Rcpp::List components(prep["components"]);
    model.component_mode = detail::as_string(components["mode"], "fixed");
    model.component_reference = detail::as_string(components["reference"]);

    const auto component_ids = detail::as_string_vector(components["ids"]);
    Rcpp::NumericVector weights = components.containsElementNamed("weights")
                                      ? Rcpp::NumericVector(components["weights"])
                                      : Rcpp::NumericVector();
    Rcpp::List attrs = components.containsElementNamed("attrs")
                           ? Rcpp::List(components["attrs"])
                           : Rcpp::List();

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
      component.weight = i < static_cast<std::size_t>(weights.size())
                             ? static_cast<double>(weights[static_cast<R_xlen_t>(i)])
                             : 1.0;
      if (component.id == "__default__") {
        component.active_leaf_indices.resize(model.leaves.size());
        for (std::size_t j = 0; j < model.leaves.size(); ++j) {
          component.active_leaf_indices[j] = static_cast<semantic::Index>(j);
        }
      } else {
        component.active_leaf_indices = active_map[component.id];
      }
      if (attrs.containsElementNamed(component.id.c_str())) {
        Rcpp::List attr(attrs[component.id]);
        if (detail::has_named_element(attr, "weight_param")) {
          component.weight_name = detail::as_string(attr["weight_param"]);
        }
        if (detail::has_named_element(attr, "n_outcomes") &&
            !Rf_isNull(attr["n_outcomes"])) {
          component.n_outcomes_override = detail::as_int(attr["n_outcomes"], 0);
        }
      }
      model.components.push_back(std::move(component));
    }
  }

  if (prep.containsElementNamed("observation") && !Rf_isNull(prep["observation"])) {
    Rcpp::List observation(prep["observation"]);
    model.observation.mode = semantic::ObservationMode::TopK;
    model.observation.n_outcomes = detail::as_int(observation["n_outcomes"], 1);
    model.observation.global_n_outcomes =
        detail::as_int(observation["global_n_outcomes"],
                       model.observation.n_outcomes);
  }

  return model;
}

} // namespace accumulatr::compile
