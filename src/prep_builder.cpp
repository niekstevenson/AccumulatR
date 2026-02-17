#include "prep_builder.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "accumulator.h"

namespace uuber {
namespace {

inline double clamp_probability(double value) {
  if (!std::isfinite(value))
    return 0.0;
  if (value < 0.0)
    return 0.0;
  if (value > 1.0)
    return 1.0;
  return value;
}

std::vector<std::string> extract_string_vector(SEXP obj) {
  if (Rf_isNull(obj))
    return {};
  Rcpp::CharacterVector chr(obj);
  std::vector<std::string> out;
  out.reserve(chr.size());
  for (R_xlen_t i = 0; i < chr.size(); ++i) {
    if (chr[i] == NA_STRING)
      continue;
    out.push_back(Rcpp::as<std::string>(chr[i]));
  }
  return out;
}

int extract_int(SEXP obj, int default_value);

std::vector<int> extract_int_vector(SEXP obj, bool subtract_one = false) {
  if (Rf_isNull(obj))
    return {};
  if (TYPEOF(obj) == VECSXP) {
    Rcpp::List lst(obj);
    std::vector<int> out;
    out.reserve(lst.size());
    for (R_xlen_t i = 0; i < lst.size(); ++i) {
      int val = extract_int(lst[i], NA_INTEGER);
      if (val == NA_INTEGER)
        continue;
      out.push_back(subtract_one ? val - 1 : val);
    }
    return out;
  }
  Rcpp::IntegerVector vec(obj);
  std::vector<int> out;
  out.reserve(vec.size());
  for (int val : vec) {
    if (val == NA_INTEGER)
      continue;
    out.push_back(subtract_one ? val - 1 : val);
  }
  return out;
}

std::vector<double> extract_double_vector(SEXP obj) {
  if (Rf_isNull(obj))
    return {};
  Rcpp::NumericVector vec(obj);
  std::vector<double> out;
  out.reserve(vec.size());
  for (double val : vec) {
    out.push_back(val);
  }
  return out;
}

std::string extract_string(SEXP obj) {
  if (Rf_isNull(obj) || obj == NA_STRING)
    return std::string();
  if (TYPEOF(obj) == VECSXP) {
    Rcpp::List lst(obj);
    if (lst.size() == 0)
      return std::string();
    return extract_string(lst[0]);
  }
  if (Rf_isString(obj)) {
    Rcpp::CharacterVector chr(obj);
    if (chr.size() == 0 || chr[0] == NA_STRING)
      return std::string();
    return Rcpp::as<std::string>(chr[0]);
  }
  try {
    return Rcpp::as<std::string>(obj);
  } catch (...) {
    return std::string();
  }
}

NodeKind parse_node_kind(const std::string &kind) {
  if (kind == "event")
    return NODE_EVENT;
  if (kind == "and")
    return NODE_AND;
  if (kind == "or")
    return NODE_OR;
  if (kind == "not")
    return NODE_NOT;
  if (kind == "guard")
    return NODE_GUARD;
  return NODE_UNKNOWN;
}

int extract_int(SEXP obj, int default_value = -1) {
  if (Rf_isNull(obj))
    return default_value;
  if (TYPEOF(obj) == VECSXP) {
    Rcpp::List lst(obj);
    if (lst.size() == 0)
      return default_value;
    return extract_int(lst[0], default_value);
  }
  int val = default_value;
  try {
    val = Rcpp::as<int>(obj);
  } catch (...) {
    return default_value;
  }
  if (val == NA_INTEGER)
    return default_value;
  return val;
}

std::vector<ProtoParamEntry> extract_params(SEXP obj) {
  if (Rf_isNull(obj))
    return {};
  Rcpp::List lst(obj);
  if (lst.size() == 0)
    return {};
  Rcpp::CharacterVector names = lst.names();
  std::vector<ProtoParamEntry> params;
  params.reserve(lst.size());
  for (R_xlen_t i = 0; i < lst.size(); ++i) {
    ProtoParamEntry entry;
    entry.name = Rcpp::as<std::string>(names[i]);
    SEXP val = lst[i];
    if (Rf_isLogical(val)) {
      Rcpp::LogicalVector lv(val);
      if (lv.size() == 1) {
        entry.tag = ParamValueTag::LogicalScalar;
        entry.logical_scalar = lv[0];
      } else {
        entry.tag = ParamValueTag::LogicalVector;
        entry.logical_values.assign(lv.begin(), lv.end());
      }
    } else {
      Rcpp::NumericVector nv(val);
      if (nv.size() == 1) {
        entry.tag = ParamValueTag::NumericScalar;
        entry.numeric_scalar = nv[0];
      } else {
        entry.tag = ParamValueTag::NumericVector;
        entry.numeric_values.assign(nv.begin(), nv.end());
      }
    }
    params.push_back(std::move(entry));
  }
  return params;
}

} // namespace

NativePrepProto build_prep_proto(const Rcpp::List &prep) {
  NativePrepProto proto;

  Rcpp::List acc_list = prep["accumulators"];
  if (!acc_list.isNULL()) {
    Rcpp::CharacterVector acc_names(acc_list.names());
    for (R_xlen_t i = 0; i < acc_list.size(); ++i) {
      if (acc_list[i] == R_NilValue)
        continue;
      Rcpp::List acc(acc_list[i]);
      ProtoAccumulator acc_proto;
      if (!acc_names.isNULL() && acc_names[i] != NA_STRING) {
        acc_proto.id = Rcpp::as<std::string>(acc_names[i]);
      }
      if (acc.containsElementNamed("dist"))
        acc_proto.dist = extract_string(acc["dist"]);
      if (acc.containsElementNamed("onset") && !Rf_isNull(acc["onset"])) {
        acc_proto.onset = Rcpp::as<double>(acc["onset"]);
      }
      acc_proto.onset_kind = ONSET_ABSOLUTE;
      acc_proto.onset_lag = 0.0;
      if (acc.containsElementNamed("onset_spec") && !Rf_isNull(acc["onset_spec"])) {
        Rcpp::List onset_spec(acc["onset_spec"]);
        std::string onset_kind = onset_spec.containsElementNamed("kind")
                                     ? extract_string(onset_spec["kind"])
                                     : std::string("absolute");
        if (onset_kind == "after") {
          std::string source = onset_spec.containsElementNamed("source")
                                   ? extract_string(onset_spec["source"])
                                   : std::string();
          std::string source_kind =
              onset_spec.containsElementNamed("source_kind")
                  ? extract_string(onset_spec["source_kind"])
                  : std::string();
          double lag = onset_spec.containsElementNamed("lag") &&
                               !Rf_isNull(onset_spec["lag"])
                           ? Rcpp::as<double>(onset_spec["lag"])
                           : 0.0;
          if (source_kind == "accumulator") {
            acc_proto.onset_kind = ONSET_AFTER_ACCUMULATOR;
          } else if (source_kind == "pool") {
            acc_proto.onset_kind = ONSET_AFTER_POOL;
          } else {
            acc_proto.onset_kind = ONSET_ABSOLUTE;
          }
          acc_proto.onset_source = source;
          acc_proto.onset_source_kind = source_kind;
          acc_proto.onset_lag = lag;
        } else {
          acc_proto.onset_kind = ONSET_ABSOLUTE;
          acc_proto.onset_source.clear();
          acc_proto.onset_source_kind.clear();
          acc_proto.onset_lag = 0.0;
        }
      }
      if (acc.containsElementNamed("q") && !Rf_isNull(acc["q"])) {
        acc_proto.q = Rcpp::as<double>(acc["q"]);
      }
      if (acc.containsElementNamed("components")) {
        acc_proto.components = extract_string_vector(acc["components"]);
      }
      if (acc.containsElementNamed("shared_trigger_id")) {
        acc_proto.shared_trigger_id = extract_string(acc["shared_trigger_id"]);
      }
      if (acc.containsElementNamed("params")) {
        acc_proto.params = extract_params(acc["params"]);
      }
      proto.accumulators.push_back(std::move(acc_proto));
    }
  }

  Rcpp::List pool_list = prep["pools"];
  if (!pool_list.isNULL()) {
    Rcpp::CharacterVector pool_names(pool_list.names());
    for (R_xlen_t i = 0; i < pool_list.size(); ++i) {
      if (pool_list[i] == R_NilValue)
        continue;
      Rcpp::List pool(pool_list[i]);
      ProtoPool pool_proto;
      if (!pool_names.isNULL() && pool_names[i] != NA_STRING) {
        pool_proto.id = Rcpp::as<std::string>(pool_names[i]);
      }
      if (pool.containsElementNamed("k") && !Rf_isNull(pool["k"])) {
        pool_proto.k = Rcpp::as<int>(pool["k"]);
      }
      if (pool.containsElementNamed("members")) {
        pool_proto.members = extract_string_vector(pool["members"]);
      }
      proto.pools.push_back(std::move(pool_proto));
    }
  }

  Rcpp::List id_index = prep[".id_index"];
  if (!id_index.isNULL()) {
    Rcpp::CharacterVector names(id_index.names());
    for (R_xlen_t i = 0; i < id_index.size(); ++i) {
      if (names.isNULL() || names[i] == NA_STRING)
        continue;
      int value = extract_int(id_index[i], NA_INTEGER);
      if (value == NA_INTEGER)
        continue;
      ProtoLabelEntry entry;
      entry.label = Rcpp::as<std::string>(names[i]);
      entry.id = value;
      proto.label_index.push_back(std::move(entry));
    }
  }

  Rcpp::List expr_compiled = prep[".expr_compiled"];
  if (!expr_compiled.isNULL()) {
    Rcpp::List node_list = expr_compiled["nodes"];
    proto.nodes.reserve(node_list.size());
    for (R_xlen_t i = 0; i < node_list.size(); ++i) {
      if (node_list[i] == R_NilValue)
        continue;
      Rcpp::List node(node_list[i]);
      ProtoNode proto_node;
      proto_node.id = extract_int(node["id"], static_cast<int>(i + 1));
      if (node.containsElementNamed("kind"))
        proto_node.kind = extract_string(node["kind"]);
      if (node.containsElementNamed("source"))
        proto_node.source = extract_string(node["source"]);
      if (node.containsElementNamed("args"))
        proto_node.args = extract_int_vector(node["args"]);
      if (node.containsElementNamed("sources"))
        proto_node.source_ids = extract_int_vector(node["sources"]);
      if (node.containsElementNamed("reference_id"))
        proto_node.reference_id = extract_int(node["reference_id"]);
      if (node.containsElementNamed("blocker_id"))
        proto_node.blocker_id = extract_int(node["blocker_id"]);

      if (node.containsElementNamed("arg"))
        proto_node.arg_id = extract_int(node["arg"]);
      if (node.containsElementNamed("needs_forced")) {
        proto_node.needs_forced = Rcpp::as<bool>(node["needs_forced"]);
      }
      if (node.containsElementNamed("scenario_sensitive")) {
        proto_node.scenario_sensitive =
            Rcpp::as<bool>(node["scenario_sensitive"]);
      }
      if (proto_node.source.empty() && node.containsElementNamed("expr")) {
        Rcpp::List expr = node["expr"];
        if (!expr.isNULL() && expr.containsElementNamed("source")) {
          proto_node.source = extract_string(expr["source"]);
        }
      }
      proto.nodes.push_back(std::move(proto_node));
    }
  }

  Rcpp::List outcome_list = prep["outcomes"];
  if (!outcome_list.isNULL()) {
    Rcpp::CharacterVector outcome_names(outcome_list.names());
    proto.outcomes.reserve(outcome_list.size());
    for (R_xlen_t i = 0; i < outcome_list.size(); ++i) {
      if (outcome_list[i] == R_NilValue)
        continue;
      if (outcome_names.isNULL() || outcome_names[i] == NA_STRING)
        continue;
      ProtoOutcome outcome_proto;
      outcome_proto.label = Rcpp::as<std::string>(outcome_names[i]);
      Rcpp::List def(outcome_list[i]);
      if (def.containsElementNamed("expr")) {
        Rcpp::RObject expr = def["expr"];
        if (!expr.isNULL()) {
          Rcpp::RObject lik_id = expr.attr(".lik_id");
          if (!lik_id.isNULL()) {
            outcome_proto.node_id = extract_int(lik_id, -1);
          }
        }
      }
      Rcpp::List options = def.containsElementNamed("options")
                               ? Rcpp::List(def["options"])
                               : Rcpp::List();
      if (!options.isNULL()) {
        if (options.containsElementNamed("component")) {
          outcome_proto.allowed_components =
              extract_string_vector(options["component"]);
        }
        if (options.containsElementNamed("map_outcome_to")) {
          Rcpp::RObject map_obj = options["map_outcome_to"];
          if (!map_obj.isNULL()) {
            Rcpp::CharacterVector map_cv(map_obj);
            if (map_cv.size() > 0) {
              Rcpp::String map_val = map_cv[0];
              if (map_val == NA_STRING) {
                outcome_proto.maps_to_na = true;
              } else {
                outcome_proto.map_target = Rcpp::as<std::string>(
                    Rcpp::CharacterVector::create(map_val)[0]);
              }
            }
          }
        }
        if (options.containsElementNamed("guess")) {
          Rcpp::RObject guess_obj = options["guess"];
          if (!guess_obj.isNULL()) {
            Rcpp::List guess_list(guess_obj);
            std::string rt_policy = "keep";
            if (guess_list.containsElementNamed("rt_policy")) {
              Rcpp::RObject pol = guess_list["rt_policy"];
              if (!pol.isNULL()) {
                rt_policy = extract_string(pol);
              }
            }
            Rcpp::CharacterVector labels =
                guess_list.containsElementNamed("labels")
                    ? Rcpp::CharacterVector(guess_list["labels"])
                    : Rcpp::CharacterVector();
            Rcpp::NumericVector weights =
                guess_list.containsElementNamed("weights")
                    ? Rcpp::NumericVector(guess_list["weights"])
                    : Rcpp::NumericVector();
            if (labels.size() == weights.size()) {
              for (R_xlen_t j = 0; j < labels.size(); ++j) {
                if (labels[j] == NA_STRING)
                  continue;
                ProtoOutcomeGuessDonor donor;
                donor.label = Rcpp::as<std::string>(labels[j]);
                donor.weight = weights[j];
                donor.rt_policy = rt_policy;
                outcome_proto.guess_donors.push_back(std::move(donor));
              }
            }
          }
        }
      }
      proto.outcomes.push_back(std::move(outcome_proto));
    }
  }

  auto parse_competitor_node_ids = [](SEXP comp_obj) -> std::vector<int> {
    if (Rf_isNull(comp_obj))
      return {};
    Rcpp::List comp_list(comp_obj);
    std::vector<int> comp_ids;
    comp_ids.reserve(comp_list.size());
    for (R_xlen_t j = 0; j < comp_list.size(); ++j) {
      Rcpp::RObject entry = comp_list[j];
      if (entry.isNULL())
        continue;
      Rcpp::RObject lik_id = entry.attr(".lik_id");
      if (lik_id.isNULL())
        continue;
      int lik_val = extract_int(lik_id, NA_INTEGER);
      if (lik_val == NA_INTEGER)
        continue;
      comp_ids.push_back(lik_val);
    }
    return comp_ids;
  };

  Rcpp::List competitors = prep[".competitors"];
  if (!competitors.isNULL() && !proto.outcomes.empty()) {
    bool by_index = false;
    Rcpp::RObject by_index_attr = competitors.attr("by_index");
    if (!by_index_attr.isNULL()) {
      by_index = Rcpp::as<bool>(by_index_attr);
    }
    if (by_index &&
        competitors.size() == static_cast<R_xlen_t>(proto.outcomes.size())) {
      for (R_xlen_t i = 0; i < competitors.size(); ++i) {
        if (Rf_isNull(competitors[i]))
          continue;
        proto.outcomes[static_cast<std::size_t>(i)].competitor_ids =
            parse_competitor_node_ids(competitors[i]);
      }
    } else {
      Rcpp::CharacterVector comp_names = competitors.names();
      for (R_xlen_t i = 0; i < competitors.size(); ++i) {
        if (Rf_isNull(competitors[i]))
          continue;
        if (comp_names.isNULL() || comp_names[i] == NA_STRING)
          continue;
        std::string outcome_label = Rcpp::as<std::string>(comp_names[i]);
        std::vector<int> comp_ids = parse_competitor_node_ids(competitors[i]);
        for (auto &outcome : proto.outcomes) {
          if (outcome.label == outcome_label) {
            outcome.competitor_ids = comp_ids;
          }
        }
      }
    }
  }

  Rcpp::List components = prep["components"];
  if (!components.isNULL()) {
    Rcpp::CharacterVector comp_ids =
        components.containsElementNamed("ids")
            ? Rcpp::CharacterVector(components["ids"])
            : Rcpp::CharacterVector();
    Rcpp::NumericVector comp_weights =
        components.containsElementNamed("weights")
            ? Rcpp::NumericVector(components["weights"])
            : Rcpp::NumericVector();
    Rcpp::List comp_attrs = components.containsElementNamed("attrs")
                                ? Rcpp::List(components["attrs"])
                                : Rcpp::List();
    proto.components.reserve(comp_ids.size());
    for (R_xlen_t i = 0; i < comp_ids.size(); ++i) {
      if (comp_ids[i] == NA_STRING)
        continue;
      ProtoComponent component_proto;
      component_proto.id = Rcpp::as<std::string>(comp_ids[i]);
      double weight = (i < comp_weights.size()) ? comp_weights[i] : 1.0;
      if (!std::isfinite(weight) || weight < 0.0)
        weight = 1.0;
      component_proto.weight = weight;
      if (!comp_attrs.isNULL() && i < comp_attrs.size() &&
          comp_attrs[i] != R_NilValue) {
        Rcpp::List attrs(comp_attrs[i]);
        if (attrs.containsElementNamed("guess")) {
          Rcpp::RObject guess_obj = attrs["guess"];
          if (!guess_obj.isNULL()) {
            Rcpp::List guess_list(guess_obj);
            component_proto.has_guess = true;
            if (guess_list.containsElementNamed("outcome")) {
              component_proto.guess.target = extract_string(guess_list["outcome"]);
            }
            if (component_proto.guess.target.empty()) {
              component_proto.guess.target = "__GUESS__";
            }
            if (guess_list.containsElementNamed("weights")) {
              Rcpp::NumericVector weights(guess_list["weights"]);
              Rcpp::CharacterVector weight_names = weights.names();
              if (!weight_names.isNULL()) {
                for (R_xlen_t j = 0; j < weights.size(); ++j) {
                  if (weight_names[j] == NA_STRING)
                    continue;
                  component_proto.guess.weight_labels.push_back(
                      Rcpp::as<std::string>(weight_names[j]));
                  component_proto.guess.weights.push_back(weights[j]);
                }
              } else {
                component_proto.guess.weights = extract_double_vector(weights);
              }
            }
          }
        }
      }
      proto.components.push_back(std::move(component_proto));
    }
  }

  return proto;
}

std::unique_ptr<NativeContext>
build_context_from_proto(const NativePrepProto &proto) {
  auto ctx = std::make_unique<NativeContext>();

  for (const auto &label : proto.label_index) {
    if (!label.label.empty()) {
      ctx->label_to_id[label.label] = label.id;
    }
  }

  ctx->accumulators.reserve(proto.accumulators.size());
  for (const auto &acc_proto : proto.accumulators) {
    NativeAccumulator acc;
    acc.id = acc_proto.id;
    acc.dist = acc_proto.dist;
    acc.onset = acc_proto.onset;
    acc.onset_kind = acc_proto.onset_kind;
    acc.onset_lag = acc_proto.onset_lag;
    acc.onset_source_acc_idx = -1;
    acc.onset_source_pool_idx = -1;
    acc.q = clamp_probability(acc_proto.q);
    acc.components = acc_proto.components;
    acc.shared_trigger_id = acc_proto.shared_trigger_id;
    try {
      acc.dist_cfg = resolve_acc_params_entries(acc.dist, acc_proto.params);
    } catch (const std::exception &e) {
      Rcpp::stop("native_context_build: accumulator '%s' %s", acc.id.c_str(),
                 e.what());
    }
    ctx->accumulators.push_back(acc);
  }

  for (std::size_t i = 0; i < ctx->accumulators.size(); ++i) {
    const auto &acc = ctx->accumulators[i];
    if (!acc.id.empty()) {
      ctx->accumulator_index[acc.id] = static_cast<int>(i);
    }
    if (!acc.shared_trigger_id.empty()) {
      ctx->shared_trigger_map[acc.shared_trigger_id].push_back(
          static_cast<int>(i));
    }
  }
  ctx->accumulator_label_ids.assign(ctx->accumulators.size(), NA_INTEGER);
  for (std::size_t i = 0; i < ctx->accumulators.size(); ++i) {
    const std::string &label = ctx->accumulators[i].id;
    auto it = ctx->label_to_id.find(label);
    if (it != ctx->label_to_id.end()) {
      ctx->accumulator_label_ids[i] = it->second;
    }
  }
  ctx->shared_trigger_groups.clear();
  ctx->shared_trigger_groups.reserve(ctx->shared_trigger_map.size());
  for (const auto &entry : ctx->shared_trigger_map) {
    SharedTriggerGroup group;
    group.acc_indices = entry.second;
    if (!group.acc_indices.empty()) {
      group.q_acc_idx = group.acc_indices.front();
    }
    ctx->shared_trigger_groups.push_back(std::move(group));
  }

  ctx->pools.reserve(proto.pools.size());
  for (const auto &pool_proto : proto.pools) {
    NativePool pool;
    pool.id = pool_proto.id;
    pool.k = pool_proto.k;
    pool.members = pool_proto.members;
    ctx->pools.push_back(pool);
  }
  for (std::size_t i = 0; i < ctx->pools.size(); ++i) {
    if (!ctx->pools[i].id.empty()) {
      ctx->pool_index[ctx->pools[i].id] = static_cast<int>(i);
    }
  }
  ctx->pool_label_ids.assign(ctx->pools.size(), NA_INTEGER);
  for (std::size_t i = 0; i < ctx->pools.size(); ++i) {
    const std::string &label = ctx->pools[i].id;
    auto it = ctx->label_to_id.find(label);
    if (it != ctx->label_to_id.end()) {
      ctx->pool_label_ids[i] = it->second;
    }
  }

  for (std::size_t i = 0; i < ctx->accumulators.size(); ++i) {
    NativeAccumulator &acc = ctx->accumulators[i];
    if (i >= proto.accumulators.size()) {
      continue;
    }
    const ProtoAccumulator &acc_proto = proto.accumulators[i];
    if (acc.onset_kind == ONSET_ABSOLUTE) {
      continue;
    }
    if (acc_proto.onset_source.empty()) {
      Rcpp::stop("native_context_build: accumulator '%s' has chained onset with empty source",
                 acc.id.c_str());
    }
    if (acc.onset_kind == ONSET_AFTER_ACCUMULATOR) {
      auto it = ctx->accumulator_index.find(acc_proto.onset_source);
      if (it == ctx->accumulator_index.end()) {
        Rcpp::stop("native_context_build: accumulator '%s' references unknown onset accumulator '%s'",
                   acc.id.c_str(), acc_proto.onset_source.c_str());
      }
      acc.onset_source_acc_idx = it->second;
      ctx->has_chained_onsets = true;
      continue;
    }
    if (acc.onset_kind == ONSET_AFTER_POOL) {
      auto it = ctx->pool_index.find(acc_proto.onset_source);
      if (it == ctx->pool_index.end()) {
        Rcpp::stop("native_context_build: accumulator '%s' references unknown onset pool '%s'",
                   acc.id.c_str(), acc_proto.onset_source.c_str());
      }
      acc.onset_source_pool_idx = it->second;
      ctx->has_chained_onsets = true;
      continue;
    }
    Rcpp::stop("native_context_build: accumulator '%s' has unsupported onset kind",
               acc.id.c_str());
  }

  auto resolve_label_ref = [&](const std::string &label) -> LabelRef {
    LabelRef ref;
    if (label.empty())
      return ref;
    auto it_label = ctx->label_to_id.find(label);
    if (it_label != ctx->label_to_id.end()) {
      ref.label_id = it_label->second;
    }
    auto it_acc = ctx->accumulator_index.find(label);
    if (it_acc != ctx->accumulator_index.end()) {
      ref.acc_idx = it_acc->second;
    }
    auto it_pool = ctx->pool_index.find(label);
    if (it_pool != ctx->pool_index.end()) {
      ref.pool_idx = it_pool->second;
    }
    // Outcome index? Not usually populated in Proto flow but handled for
    // completeness if it were? Actually proto doesn't have outcomes. But
    // assuming usage consistency:
    auto it_out = ctx->outcome_index.find(label);
    if (it_out != ctx->outcome_index.end() && !it_out->second.empty()) {
      ref.outcome_idx = it_out->second[0];
    }
    return ref;
  };

  for (auto &pool : ctx->pools) {
    pool.member_refs.clear();
    pool.member_refs.reserve(pool.members.size());
    for (const auto &member : pool.members) {
      pool.member_refs.push_back(resolve_label_ref(member));
    }
  }

  ctx->nodes.reserve(proto.nodes.size());
  for (const auto &node_proto : proto.nodes) {
    NativeNode node;
    node.id = node_proto.id;
    node.kind = node_proto.kind;
    node.kind_id = parse_node_kind(node_proto.kind);
    node.source = node_proto.source;
    node.source_ref = resolve_label_ref(node.source);
    node.args = node_proto.args;
    node.source_ids = node_proto.source_ids;
    if (!node.source_ids.empty()) {
      std::sort(node.source_ids.begin(), node.source_ids.end());
      node.source_ids.erase(
          std::unique(node.source_ids.begin(), node.source_ids.end()),
          node.source_ids.end());
    }
    node.reference_id = node_proto.reference_id;
    node.blocker_id = node_proto.blocker_id;

    node.arg_id = node_proto.arg_id;
    node.needs_forced = node_proto.needs_forced;
    node.scenario_sensitive = node_proto.scenario_sensitive;
    ctx->node_index[node.id] = static_cast<int>(ctx->nodes.size());
    ctx->nodes.push_back(std::move(node));
  }

  ctx->outcome_index.clear();
  ctx->outcome_labels.clear();
  ctx->outcome_label_ids.clear();
  ctx->outcome_info.clear();
  ctx->outcome_labels.reserve(proto.outcomes.size());
  ctx->outcome_label_ids.reserve(proto.outcomes.size());
  ctx->outcome_info.reserve(proto.outcomes.size());
  for (const auto &outcome_proto : proto.outcomes) {
    int idx = static_cast<int>(ctx->outcome_info.size());
    ctx->outcome_index[outcome_proto.label].push_back(idx);
    ctx->outcome_labels.push_back(outcome_proto.label);
    auto it_label = ctx->label_to_id.find(outcome_proto.label);
    ctx->outcome_label_ids.push_back(
        (it_label == ctx->label_to_id.end()) ? NA_INTEGER : it_label->second);
    OutcomeContextInfo info;
    info.node_id = outcome_proto.node_id;
    info.competitor_ids = outcome_proto.competitor_ids;
    info.allowed_components = outcome_proto.allowed_components;
    info.maps_to_na = outcome_proto.maps_to_na;
    ctx->outcome_info.push_back(std::move(info));
  }
  for (std::size_t i = 0; i < proto.outcomes.size(); ++i) {
    const ProtoOutcome &outcome_proto = proto.outcomes[i];
    int source_idx = static_cast<int>(i);
    if (source_idx >= static_cast<int>(ctx->outcome_info.size()))
      break;
    if (!outcome_proto.map_target.empty()) {
      auto it = ctx->outcome_index.find(outcome_proto.map_target);
      if (it != ctx->outcome_index.end()) {
        for (int target_idx : it->second) {
          if (target_idx >= 0 &&
              target_idx < static_cast<int>(ctx->outcome_info.size())) {
            ctx->outcome_info[static_cast<std::size_t>(target_idx)]
                .alias_sources.push_back(source_idx);
          }
        }
      }
    }
    for (const auto &guess_proto : outcome_proto.guess_donors) {
      auto it = ctx->outcome_index.find(guess_proto.label);
      if (it == ctx->outcome_index.end())
        continue;
      for (int target_idx : it->second) {
        if (target_idx < 0 ||
            target_idx >= static_cast<int>(ctx->outcome_info.size())) {
          continue;
        }
        OutcomeGuessDonor donor;
        donor.outcome_idx = source_idx;
        donor.weight = guess_proto.weight;
        donor.rt_policy =
            guess_proto.rt_policy.empty() ? "keep" : guess_proto.rt_policy;
        ctx->outcome_info[static_cast<std::size_t>(target_idx)]
            .guess_donors.push_back(std::move(donor));
      }
    }
  }

  ctx->component_index.clear();
  ctx->component_info.clear();
  ctx->components.ids.clear();
  ctx->components.base_weights.clear();
  ctx->components.leader_idx.clear();
  ctx->components.ids.reserve(proto.components.size());
  ctx->components.base_weights.reserve(proto.components.size());
  for (const auto &component_proto : proto.components) {
    int idx = static_cast<int>(ctx->component_info.size());
    ctx->component_index[component_proto.id] = idx;
    ComponentContextInfo meta;
    if (component_proto.has_guess) {
      meta.guess.target =
          component_proto.guess.target.empty() ? "__GUESS__"
                                               : component_proto.guess.target;
      if (!meta.guess.target.empty() && meta.guess.target != "__GUESS__") {
        auto it = ctx->outcome_index.find(meta.guess.target);
        if (it != ctx->outcome_index.end() && !it->second.empty()) {
          meta.guess.target_outcome_idx = it->second[0];
        }
      }
      std::size_t n =
          std::min(component_proto.guess.weight_labels.size(),
                   component_proto.guess.weights.size());
      for (std::size_t wi = 0; wi < n; ++wi) {
        const std::string &label = component_proto.guess.weight_labels[wi];
        const double keep_prob = component_proto.guess.weights[wi];
        if (label == "<NA>" || label == "NA") {
          meta.guess.keep_weight_na = keep_prob;
          meta.guess.has_keep_weight_na = true;
          continue;
        }
        auto it = ctx->outcome_index.find(label);
        if (it == ctx->outcome_index.end())
          continue;
        for (int out_idx : it->second) {
          meta.guess.keep_weights.emplace_back(out_idx, keep_prob);
        }
      }
      meta.guess.valid =
          !meta.guess.keep_weights.empty() || meta.guess.has_keep_weight_na;
    }
    ctx->component_info.push_back(std::move(meta));
    ctx->components.ids.push_back(component_proto.id);
    double weight = component_proto.weight;
    if (!std::isfinite(weight) || weight < 0.0)
      weight = 1.0;
    ctx->components.base_weights.push_back(weight);
  }
  if (ctx->components.ids.empty()) {
    ctx->components.ids.push_back("__default__");
    ctx->components.base_weights.push_back(1.0);
  }
  ctx->components.leader_idx.assign(ctx->components.ids.size(), -1);
  for (std::size_t c = 0; c < ctx->components.ids.size(); ++c) {
    const std::string &cid = ctx->components.ids[c];
    for (std::size_t a = 0; a < ctx->accumulators.size(); ++a) {
      const auto &comps = ctx->accumulators[a].components;
      if (std::find(comps.begin(), comps.end(), cid) != comps.end()) {
        ctx->components.leader_idx[c] = static_cast<int>(a);
        break;
      }
    }
    if (ctx->components.leader_idx[c] < 0 && !ctx->accumulators.empty()) {
      ctx->components.leader_idx[c] = 0;
    }
  }

  for (auto &acc : ctx->accumulators) {
    acc.component_indices.clear();
    acc.component_indices.reserve(acc.components.size());
    for (const auto &comp : acc.components) {
      auto it = ctx->component_index.find(comp);
      if (it != ctx->component_index.end()) {
        acc.component_indices.push_back(it->second);
      }
    }
  }
  for (auto &node : ctx->nodes) {
    if (node.source.empty())
      continue;
    auto it = ctx->outcome_index.find(node.source);
    if (it != ctx->outcome_index.end() && !it->second.empty()) {
      node.source_ref.outcome_idx = it->second[0];
    }
  }
  for (auto &pool : ctx->pools) {
    if (pool.members.empty() || pool.member_refs.empty())
      continue;
    std::size_t n = pool.members.size();
    if (pool.member_refs.size() != n) {
      pool.member_refs.resize(n);
    }
    for (std::size_t i = 0; i < n; ++i) {
      auto it = ctx->outcome_index.find(pool.members[i]);
      if (it != ctx->outcome_index.end() && !it->second.empty()) {
        pool.member_refs[i].outcome_idx = it->second[0];
      }
    }
  }

  build_ir_context(*ctx);

  return ctx;
}

} // namespace uuber
