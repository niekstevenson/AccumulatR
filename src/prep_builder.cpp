#include "prep_builder.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#include "accumulator.h"
#include "kernel_program.h"

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

void validate_guard_transition_completeness(const NativeContext &ctx) {
  const IrContext &ir = ctx.ir;
  const KernelStateGraph &graph = ctx.kernel_state_graph;
  if (graph.node_guard_transition_idx.size() != ir.nodes.size()) {
    Rcpp::stop("IR guard transition index map size mismatch");
  }
  for (int node_idx = 0; node_idx < static_cast<int>(ir.nodes.size());
       ++node_idx) {
    const IrNode &node = ir.nodes[static_cast<std::size_t>(node_idx)];
    if (node.source_id_count <= 0) {
      continue;
    }
    const int tr_idx =
        graph.node_guard_transition_idx[static_cast<std::size_t>(node_idx)];
    if (tr_idx < 0 ||
        tr_idx >= static_cast<int>(graph.guard_transitions.size())) {
      Rcpp::stop(
          "IR guard transition missing for node_id=%d node_idx=%d source_count=%d",
          node.node_id, node_idx, node.source_id_count);
    }
    const KernelGuardTransition &tr =
        graph.guard_transitions[static_cast<std::size_t>(tr_idx)];
    if (tr.node_idx != node_idx) {
      Rcpp::stop(
          "IR guard transition node mismatch for node_id=%d node_idx=%d source_count=%d",
          node.node_id, node_idx, node.source_id_count);
    }
    if (tr.source_mask_begin < 0 || tr.source_mask_count <= 0 ||
        tr.source_mask_begin + tr.source_mask_count >
            static_cast<int>(ir.node_source_masks.size())) {
      Rcpp::stop(
          "IR guard transition mask invalid for node_id=%d node_idx=%d source_count=%d",
          node.node_id, node_idx, node.source_id_count);
    }
  }
}

IrNodeOp parse_node_op(const std::string &kind, int pool_idx) {
  if (kind == "event") {
    return (pool_idx >= 0) ? IrNodeOp::EventPool : IrNodeOp::EventAcc;
  }
  if (kind == "and")
    return IrNodeOp::And;
  if (kind == "or")
    return IrNodeOp::Or;
  if (kind == "not")
    return IrNodeOp::Not;
  if (kind == "guard")
    return IrNodeOp::Guard;
  return IrNodeOp::And;
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

  std::unordered_map<std::string, int> label_id_by_name;
  label_id_by_name.reserve(proto.label_index.size());
  for (const auto &entry : proto.label_index) {
    label_id_by_name[entry.label] = entry.id;
  }
  std::unordered_map<std::string, int> acc_idx_by_name;
  acc_idx_by_name.reserve(proto.accumulators.size());
  for (std::size_t i = 0; i < proto.accumulators.size(); ++i) {
    if (!proto.accumulators[i].id.empty()) {
      acc_idx_by_name[proto.accumulators[i].id] = static_cast<int>(i);
    }
  }
  std::unordered_map<std::string, int> pool_idx_by_name;
  pool_idx_by_name.reserve(proto.pools.size());
  for (std::size_t i = 0; i < proto.pools.size(); ++i) {
    if (!proto.pools[i].id.empty()) {
      pool_idx_by_name[proto.pools[i].id] = static_cast<int>(i);
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
      std::string kind =
          node.containsElementNamed("kind") ? extract_string(node["kind"])
                                            : std::string();
      std::string source =
          node.containsElementNamed("source") ? extract_string(node["source"])
                                              : std::string();
      if (source.empty() && node.containsElementNamed("expr")) {
        Rcpp::List expr = node["expr"];
        if (!expr.isNULL() && expr.containsElementNamed("source")) {
          source = extract_string(expr["source"]);
        }
      }
      std::vector<int> child_ids;
      if (kind == "not") {
        if (node.containsElementNamed("arg")) {
          int arg_id = extract_int(node["arg"]);
          if (arg_id >= 0) {
            child_ids.push_back(arg_id);
          }
        }
      } else if (node.containsElementNamed("args")) {
        child_ids = extract_int_vector(node["args"]);
      }
      if (!child_ids.empty()) {
        std::sort(child_ids.begin(), child_ids.end());
        child_ids.erase(std::unique(child_ids.begin(), child_ids.end()),
                        child_ids.end());
      }
      proto_node.children = std::move(child_ids);
      if (node.containsElementNamed("sources")) {
        proto_node.source_label_ids = extract_int_vector(node["sources"]);
        std::sort(proto_node.source_label_ids.begin(),
                  proto_node.source_label_ids.end());
        proto_node.source_label_ids.erase(
            std::unique(proto_node.source_label_ids.begin(),
                        proto_node.source_label_ids.end()),
            proto_node.source_label_ids.end());
      }
      if (node.containsElementNamed("reference_id"))
        proto_node.reference_id = extract_int(node["reference_id"]);
      if (node.containsElementNamed("blocker_id"))
        proto_node.blocker_id = extract_int(node["blocker_id"]);

      std::uint32_t flags = IR_NODE_FLAG_NONE;
      if (node.containsElementNamed("needs_forced")) {
        if (Rcpp::as<bool>(node["needs_forced"])) {
          flags |= IR_NODE_FLAG_NEEDS_FORCED;
        }
      }
      if (node.containsElementNamed("scenario_sensitive")) {
        if (Rcpp::as<bool>(node["scenario_sensitive"])) {
          flags |= IR_NODE_FLAG_SCENARIO_SENSITIVE;
        }
      }
      if (source == "__DEADLINE__") {
        flags |= IR_NODE_FLAG_SPECIAL_DEADLINE;
      } else if (source == "__GUESS__") {
        flags |= IR_NODE_FLAG_SPECIAL_GUESS;
      }
      proto_node.flags = flags;

      auto it_acc = acc_idx_by_name.find(source);
      if (it_acc != acc_idx_by_name.end()) {
        proto_node.event_acc_idx = it_acc->second;
      }
      auto it_pool = pool_idx_by_name.find(source);
      if (it_pool != pool_idx_by_name.end()) {
        proto_node.event_pool_idx = it_pool->second;
      }
      auto it_label = label_id_by_name.find(source);
      if (it_label != label_id_by_name.end()) {
        proto_node.event_label_id = it_label->second;
      }
      proto_node.op =
          static_cast<std::uint8_t>(parse_node_op(kind, proto_node.event_pool_idx));

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

  if (!proto.outcomes.empty()) {
    std::unordered_map<int, int> label_id_to_outcome_idx;
    label_id_to_outcome_idx.reserve(proto.outcomes.size());
    for (std::size_t oi = 0; oi < proto.outcomes.size(); ++oi) {
      auto it = label_id_by_name.find(proto.outcomes[oi].label);
      if (it == label_id_by_name.end()) {
        continue;
      }
      if (label_id_to_outcome_idx.find(it->second) ==
          label_id_to_outcome_idx.end()) {
        label_id_to_outcome_idx[it->second] = static_cast<int>(oi);
      }
    }
    for (auto &node : proto.nodes) {
      if (node.event_label_id < 0) {
        continue;
      }
      auto it = label_id_to_outcome_idx.find(node.event_label_id);
      if (it != label_id_to_outcome_idx.end()) {
        node.event_outcome_idx = it->second;
      }
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
  auto sort_unique_ids = [](std::vector<int> &ids) {
    if (ids.empty()) {
      return;
    }
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
  };

  auto mask_word_count = [](int bit_count) -> int {
    if (bit_count <= 0)
      return 0;
    return (bit_count + 63) / 64;
  };

  auto set_mask_bit = [](std::vector<std::uint64_t> &mask, int bit_idx) {
    if (bit_idx < 0)
      return;
    int word = bit_idx / 64;
    int bit = bit_idx % 64;
    if (word < 0 || word >= static_cast<int>(mask.size()))
      return;
    mask[static_cast<std::size_t>(word)] |= (1ULL << bit);
  };

  auto append_mask_words =
      [](std::vector<std::uint64_t> &flat,
         const std::vector<std::uint64_t> &mask) -> int {
    if (mask.empty()) {
      return -1;
    }
    int begin = static_cast<int>(flat.size());
    flat.insert(flat.end(), mask.begin(), mask.end());
    return begin;
  };

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
    ctx->accumulators.push_back(std::move(acc));
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
    auto it = ctx->label_to_id.find(ctx->accumulators[i].id);
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
    ctx->pools.push_back(std::move(pool));
  }
  for (std::size_t i = 0; i < ctx->pools.size(); ++i) {
    if (!ctx->pools[i].id.empty()) {
      ctx->pool_index[ctx->pools[i].id] = static_cast<int>(i);
    }
  }
  ctx->pool_label_ids.assign(ctx->pools.size(), NA_INTEGER);
  for (std::size_t i = 0; i < ctx->pools.size(); ++i) {
    auto it = ctx->label_to_id.find(ctx->pools[i].id);
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
      Rcpp::stop("native_context_build: accumulator '%s' has chained onset with "
                 "empty source",
                 acc.id.c_str());
    }
    if (acc.onset_kind == ONSET_AFTER_ACCUMULATOR) {
      auto it = ctx->accumulator_index.find(acc_proto.onset_source);
      if (it == ctx->accumulator_index.end()) {
        Rcpp::stop("native_context_build: accumulator '%s' references unknown "
                   "onset accumulator '%s'",
                   acc.id.c_str(), acc_proto.onset_source.c_str());
      }
      acc.onset_source_acc_idx = it->second;
      ctx->has_chained_onsets = true;
      continue;
    }
    if (acc.onset_kind == ONSET_AFTER_POOL) {
      auto it = ctx->pool_index.find(acc_proto.onset_source);
      if (it == ctx->pool_index.end()) {
        Rcpp::stop(
            "native_context_build: accumulator '%s' references unknown onset "
            "pool '%s'",
            acc.id.c_str(), acc_proto.onset_source.c_str());
      }
      acc.onset_source_pool_idx = it->second;
      ctx->has_chained_onsets = true;
      continue;
    }
    Rcpp::stop("native_context_build: accumulator '%s' has unsupported onset "
               "kind",
               acc.id.c_str());
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
    if (source_idx >= static_cast<int>(ctx->outcome_info.size())) {
      break;
    }
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
      if (it == ctx->outcome_index.end()) {
        continue;
      }
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
      meta.guess.target_is_guess = (meta.guess.target == "__GUESS__");
      if (!meta.guess.target.empty() && !meta.guess.target_is_guess) {
        auto it = ctx->outcome_index.find(meta.guess.target);
        if (it != ctx->outcome_index.end() && !it->second.empty()) {
          meta.guess.target_outcome_idx = it->second[0];
        }
        auto it_label = ctx->label_to_id.find(meta.guess.target);
        if (it_label != ctx->label_to_id.end()) {
          meta.guess.target_label_id = it_label->second;
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
        if (it == ctx->outcome_index.end()) {
          continue;
        }
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
    if (!std::isfinite(weight) || weight < 0.0) {
      weight = 1.0;
    }
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

  auto resolve_label_ref = [&](const std::string &label) -> LabelRef {
    LabelRef ref;
    if (label.empty()) {
      return ref;
    }
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

  IrContext ir;
  ir.valid = false;
  ir.id_to_node_idx.reserve(proto.nodes.size());
  for (std::size_t i = 0; i < proto.nodes.size(); ++i) {
    ir.id_to_node_idx[proto.nodes[i].id] = static_cast<int>(i);
  }

  std::unordered_set<int> label_ids_set;
  label_ids_set.reserve(ctx->label_to_id.size() + proto.nodes.size() * 2);
  for (const auto &kv : ctx->label_to_id) {
    if (kv.second >= 0) {
      label_ids_set.insert(kv.second);
    }
  }
  for (const auto &node : proto.nodes) {
    if (node.event_label_id >= 0) {
      label_ids_set.insert(node.event_label_id);
    }
    for (int source_id : node.source_label_ids) {
      if (source_id >= 0) {
        label_ids_set.insert(source_id);
      }
    }
  }
  std::vector<int> label_ids(label_ids_set.begin(), label_ids_set.end());
  std::sort(label_ids.begin(), label_ids.end());
  ir.label_id_to_bit_idx.reserve(label_ids.size());
  for (std::size_t i = 0; i < label_ids.size(); ++i) {
    ir.label_id_to_bit_idx[label_ids[i]] = static_cast<int>(i);
  }
  ir.source_mask_words = mask_word_count(static_cast<int>(label_ids.size()));

  const int n_components = static_cast<int>(ctx->components.ids.size());
  ir.component_mask_words = mask_word_count(n_components);
  std::vector<std::uint64_t> all_component_mask(
      static_cast<std::size_t>(ir.component_mask_words), 0ULL);
  for (int i = 0; i < n_components; ++i) {
    set_mask_bit(all_component_mask, i);
  }

  std::unordered_map<int, int> label_id_to_first_outcome;
  for (std::size_t oi = 0; oi < ctx->outcome_label_ids.size(); ++oi) {
    int label_id = ctx->outcome_label_ids[oi];
    if (label_id < 0 || label_id == NA_INTEGER) {
      continue;
    }
    if (label_id_to_first_outcome.find(label_id) ==
        label_id_to_first_outcome.end()) {
      label_id_to_first_outcome[label_id] = static_cast<int>(oi);
    }
  }

  ir.nodes.reserve(proto.nodes.size());
  for (std::size_t i = 0; i < proto.nodes.size(); ++i) {
    const ProtoNode &node = proto.nodes[i];
    IrNode ir_node;
    ir_node.op = static_cast<IrNodeOp>(node.op);
    ir_node.reference_idx = -1;
    ir_node.blocker_idx = -1;
    ir_node.flags = node.flags;
    ir_node.event_idx = -1;
    ir_node.node_id = node.id;

    if (node.reference_id >= 0) {
      auto it = ir.id_to_node_idx.find(node.reference_id);
      if (it != ir.id_to_node_idx.end()) {
        ir_node.reference_idx = it->second;
      }
    }
    if (node.blocker_id >= 0) {
      auto it = ir.id_to_node_idx.find(node.blocker_id);
      if (it != ir.id_to_node_idx.end()) {
        ir_node.blocker_idx = it->second;
      }
    }

    std::vector<int> dense_children;
    dense_children.reserve(node.children.size());
    for (int child_id : node.children) {
      auto it = ir.id_to_node_idx.find(child_id);
      if (it != ir.id_to_node_idx.end()) {
        dense_children.push_back(it->second);
      }
    }
    ir_node.child_begin =
        dense_children.empty() ? -1 : static_cast<int>(ir.node_children.size());
    ir_node.child_count = static_cast<int>(dense_children.size());
    ir.node_children.insert(ir.node_children.end(), dense_children.begin(),
                            dense_children.end());

    std::vector<int> source_ids = node.source_label_ids;
    if (node.event_label_id >= 0) {
      source_ids.push_back(node.event_label_id);
    }
    sort_unique_ids(source_ids);
    ir_node.source_id_begin = source_ids.empty()
                                  ? -1
                                  : static_cast<int>(ir.node_source_label_ids.size());
    ir_node.source_id_count = static_cast<int>(source_ids.size());
    ir.node_source_label_ids.insert(ir.node_source_label_ids.end(),
                                    source_ids.begin(), source_ids.end());

    std::vector<std::uint64_t> source_mask(
        static_cast<std::size_t>(ir.source_mask_words), 0ULL);
    if (ir.source_mask_words > 0) {
      for (int source_id : source_ids) {
        auto it = ir.label_id_to_bit_idx.find(source_id);
        if (it != ir.label_id_to_bit_idx.end()) {
          set_mask_bit(source_mask, it->second);
        }
      }
    }
    ir_node.source_mask_begin =
        append_mask_words(ir.node_source_masks, source_mask);
    ir_node.source_mask_count = ir.source_mask_words;

    if (node.op == static_cast<std::uint8_t>(IrNodeOp::EventAcc) ||
        node.op == static_cast<std::uint8_t>(IrNodeOp::EventPool)) {
      IrEvent event;
      event.node_idx = static_cast<int>(i);
      event.acc_idx = node.event_acc_idx;
      event.pool_idx = node.event_pool_idx;
      event.label_id = node.event_label_id;
      event.outcome_idx = node.event_outcome_idx;
      if (event.outcome_idx < 0 && event.label_id >= 0) {
        auto it = label_id_to_first_outcome.find(event.label_id);
        if (it != label_id_to_first_outcome.end()) {
          event.outcome_idx = it->second;
        }
      }
      std::vector<std::uint64_t> component_mask(
          static_cast<std::size_t>(ir.component_mask_words), 0ULL);
      if (ir.component_mask_words > 0) {
        if (event.acc_idx >= 0 &&
            event.acc_idx < static_cast<int>(ctx->accumulators.size()) &&
            !ctx->accumulators[static_cast<std::size_t>(event.acc_idx)]
                 .component_indices.empty()) {
          const auto &indices =
              ctx->accumulators[static_cast<std::size_t>(event.acc_idx)]
                  .component_indices;
          for (int comp_idx : indices) {
            set_mask_bit(component_mask, comp_idx);
          }
        } else {
          component_mask = all_component_mask;
        }
      }
      event.component_mask_offset =
          append_mask_words(ir.component_masks, component_mask);
      ir_node.event_idx = static_cast<int>(ir.events.size());
      ir.events.push_back(std::move(event));
    }

    ir.nodes.push_back(std::move(ir_node));
  }

  ir.outcomes.reserve(ctx->outcome_info.size());
  for (std::size_t i = 0; i < ctx->outcome_info.size(); ++i) {
    const OutcomeContextInfo &info = ctx->outcome_info[i];
    IrOutcome ir_out;
    ir_out.label_id =
        (i < ctx->outcome_label_ids.size()) ? ctx->outcome_label_ids[i] : -1;
    if (ir_out.label_id == NA_INTEGER) {
      ir_out.label_id = -1;
    }
    auto node_it = ir.id_to_node_idx.find(info.node_id);
    ir_out.node_idx = (node_it == ir.id_to_node_idx.end()) ? -1 : node_it->second;

    ir_out.competitor_begin =
        info.competitor_ids.empty() ? -1 : static_cast<int>(ir.outcome_competitors.size());
    for (int comp_id : info.competitor_ids) {
      auto it = ir.id_to_node_idx.find(comp_id);
      if (it != ir.id_to_node_idx.end()) {
        ir.outcome_competitors.push_back(it->second);
      }
    }
    ir_out.competitor_count =
        (ir_out.competitor_begin < 0)
            ? 0
            : static_cast<int>(ir.outcome_competitors.size()) -
                  ir_out.competitor_begin;

    std::vector<std::uint64_t> allowed_mask(
        static_cast<std::size_t>(ir.component_mask_words), 0ULL);
    if (ir.component_mask_words > 0) {
      if (info.allowed_components.empty()) {
        allowed_mask = all_component_mask;
      } else {
        for (const auto &component : info.allowed_components) {
          auto it = ctx->component_index.find(component);
          if (it != ctx->component_index.end()) {
            set_mask_bit(allowed_mask, it->second);
          }
        }
      }
    }
    ir_out.allowed_component_mask_offset =
        append_mask_words(ir.component_masks, allowed_mask);
    ir_out.maps_to_na = info.maps_to_na;

    ir_out.alias_begin = info.alias_sources.empty()
                             ? -1
                             : static_cast<int>(ir.outcome_alias_sources.size());
    ir.outcome_alias_sources.insert(ir.outcome_alias_sources.end(),
                                    info.alias_sources.begin(),
                                    info.alias_sources.end());
    ir_out.alias_count =
        (ir_out.alias_begin < 0)
            ? 0
            : static_cast<int>(ir.outcome_alias_sources.size()) -
                  ir_out.alias_begin;

    ir_out.guess_begin = info.guess_donors.empty()
                             ? -1
                             : static_cast<int>(ir.outcome_guess_donors.size());
    for (const auto &donor : info.guess_donors) {
      IrOutcomeGuessDonor ir_donor;
      ir_donor.outcome_idx = donor.outcome_idx;
      ir_donor.weight = donor.weight;
      if (donor.rt_policy == "keep") {
        ir_donor.rt_policy_code = 0;
      } else if (donor.rt_policy == "na") {
        ir_donor.rt_policy_code = 1;
      } else if (donor.rt_policy == "drop") {
        ir_donor.rt_policy_code = 2;
      } else {
        ir_donor.rt_policy_code = 3;
      }
      ir.outcome_guess_donors.push_back(std::move(ir_donor));
    }
    ir_out.guess_count =
        (ir_out.guess_begin < 0)
            ? 0
            : static_cast<int>(ir.outcome_guess_donors.size()) -
                  ir_out.guess_begin;

    ir.outcomes.push_back(std::move(ir_out));
    if (ir.outcomes.back().label_id >= 0) {
      ir.label_id_to_outcomes[ir.outcomes.back().label_id].push_back(
          static_cast<int>(i));
    }
    if (ir.outcomes.back().node_idx >= 0) {
      ir.node_idx_to_outcomes[ir.outcomes.back().node_idx].push_back(
          static_cast<int>(i));
    }
  }

  auto mix_lookup_hash = [](std::uint64_t hash, int value) -> std::uint64_t {
    std::uint32_t v = static_cast<std::uint32_t>(value);
    for (int i = 0; i < 4; ++i) {
      std::uint8_t b = static_cast<std::uint8_t>((v >> (8 * i)) & 0xFFU);
      hash ^= static_cast<std::uint64_t>(b);
      hash *= 1099511628211ULL;
    }
    return hash;
  };

  auto shared_gate_lookup_key = [&](int node_idx,
                                    const std::vector<int> &competitors)
      -> std::uint64_t {
    std::uint64_t hash = 1469598103934665603ULL;
    hash = mix_lookup_hash(hash, node_idx);
    hash = mix_lookup_hash(hash, static_cast<int>(competitors.size()));
    for (int comp : competitors) {
      hash = mix_lookup_hash(hash, comp);
    }
    return hash;
  };

  auto mask_overlap = [&](const IrNode &a, const IrNode &b) -> bool {
    if (a.source_mask_count <= 0 || b.source_mask_count <= 0) {
      return false;
    }
    if (a.source_mask_count != b.source_mask_count) {
      return false;
    }
    for (int i = 0; i < a.source_mask_count; ++i) {
      std::uint64_t wa = ir.node_source_masks[static_cast<std::size_t>(
          a.source_mask_begin + i)];
      std::uint64_t wb = ir.node_source_masks[static_cast<std::size_t>(
          b.source_mask_begin + i)];
      if ((wa & wb) != 0ULL) {
        return true;
      }
    }
    return false;
  };

  auto and_event_pair = [&](int node_idx, int &event_a, int &event_b) -> bool {
    if (node_idx < 0 || node_idx >= static_cast<int>(ir.nodes.size())) {
      return false;
    }
    const IrNode &node = ir.nodes[static_cast<std::size_t>(node_idx)];
    if (node.op != IrNodeOp::And || node.child_count != 2) {
      return false;
    }
    if (node.child_begin < 0 ||
        node.child_begin + 1 >= static_cast<int>(ir.node_children.size())) {
      return false;
    }
    int c0 = ir.node_children[static_cast<std::size_t>(node.child_begin)];
    int c1 = ir.node_children[static_cast<std::size_t>(node.child_begin + 1)];
    if (c0 < 0 || c1 < 0 || c0 >= static_cast<int>(ir.nodes.size()) ||
        c1 >= static_cast<int>(ir.nodes.size())) {
      return false;
    }
    const IrNode &n0 = ir.nodes[static_cast<std::size_t>(c0)];
    const IrNode &n1 = ir.nodes[static_cast<std::size_t>(c1)];
    if (n0.event_idx < 0 || n1.event_idx < 0) {
      return false;
    }
    event_a = n0.event_idx;
    event_b = n1.event_idx;
    return true;
  };

  auto detect_shared_gate_pair_ir = [&](int node_idx, int competitor_idx,
                                        IrSharedGateSpec &spec_out) -> bool {
    int ta = -1, tb = -1, ca = -1, cb = -1;
    if (!and_event_pair(node_idx, ta, tb) ||
        !and_event_pair(competitor_idx, ca, cb)) {
      return false;
    }
    if (ta < 0 || tb < 0 || ca < 0 || cb < 0) {
      return false;
    }
    if (ta >= static_cast<int>(ir.events.size()) ||
        tb >= static_cast<int>(ir.events.size()) ||
        ca >= static_cast<int>(ir.events.size()) ||
        cb >= static_cast<int>(ir.events.size())) {
      return false;
    }

    const int tl0 = ir.events[static_cast<std::size_t>(ta)].label_id;
    const int tl1 = ir.events[static_cast<std::size_t>(tb)].label_id;
    const int cl0 = ir.events[static_cast<std::size_t>(ca)].label_id;
    const int cl1 = ir.events[static_cast<std::size_t>(cb)].label_id;
    if (tl0 < 0 || tl1 < 0 || cl0 < 0 || cl1 < 0) {
      return false;
    }

    int pair_x = -1, pair_y = -1, pair_c = -1;
    if (tl0 == cl0) {
      pair_c = ta;
      pair_x = tb;
      pair_y = cb;
    } else if (tl0 == cl1) {
      pair_c = ta;
      pair_x = tb;
      pair_y = ca;
    } else if (tl1 == cl0) {
      pair_c = tb;
      pair_x = ta;
      pair_y = cb;
    } else if (tl1 == cl1) {
      pair_c = tb;
      pair_x = ta;
      pair_y = ca;
    } else {
      return false;
    }
    if (pair_x < 0 || pair_y < 0 || pair_c < 0) {
      return false;
    }
    const int x_label = ir.events[static_cast<std::size_t>(pair_x)].label_id;
    const int y_label = ir.events[static_cast<std::size_t>(pair_y)].label_id;
    const int c_label = ir.events[static_cast<std::size_t>(pair_c)].label_id;
    if (x_label < 0 || y_label < 0 || c_label < 0) {
      return false;
    }
    if (x_label == y_label || x_label == c_label || y_label == c_label) {
      return false;
    }

    spec_out = IrSharedGateSpec();
    spec_out.kind = IrSharedGateKind::Pair;
    spec_out.node_idx = node_idx;
    spec_out.pair_x_event_idx = pair_x;
    spec_out.pair_y_event_idx = pair_y;
    spec_out.pair_c_event_idx = pair_c;
    return true;
  };

  auto detect_shared_gate_nway_ir = [&](int node_idx,
                                        const std::vector<int> &competitor_nodes,
                                        std::vector<int> &nway_events,
                                        IrSharedGateSpec &spec_out) -> bool {
    if (competitor_nodes.size() < 2) {
      return false;
    }
    std::vector<int> nodes;
    nodes.reserve(1 + competitor_nodes.size());
    nodes.push_back(node_idx);
    nodes.insert(nodes.end(), competitor_nodes.begin(), competitor_nodes.end());

    struct Pair {
      int e0{-1};
      int e1{-1};
      int l0{-1};
      int l1{-1};
    };
    std::vector<Pair> pairs(nodes.size());
    for (std::size_t i = 0; i < nodes.size(); ++i) {
      int e0 = -1, e1 = -1;
      if (!and_event_pair(nodes[i], e0, e1)) {
        return false;
      }
      if (e0 < 0 || e1 < 0 || e0 >= static_cast<int>(ir.events.size()) ||
          e1 >= static_cast<int>(ir.events.size())) {
        return false;
      }
      pairs[i] =
          Pair{e0, e1, ir.events[static_cast<std::size_t>(e0)].label_id,
               ir.events[static_cast<std::size_t>(e1)].label_id};
      if (pairs[i].l0 < 0 || pairs[i].l1 < 0 || pairs[i].l0 == pairs[i].l1) {
        return false;
      }
    }

    const int cand1 = pairs[0].l0;
    const int cand2 = pairs[0].l1;
    bool cand1_all = true;
    bool cand2_all = true;
    for (std::size_t i = 1; i < pairs.size(); ++i) {
      const Pair &p = pairs[i];
      if (!(p.l0 == cand1 || p.l1 == cand1)) {
        cand1_all = false;
      }
      if (!(p.l0 == cand2 || p.l1 == cand2)) {
        cand2_all = false;
      }
    }
    if (cand1_all == cand2_all) {
      return false;
    }
    const int gate_label = cand1_all ? cand1 : cand2;

    int gate_event_idx = -1;
    int target_event_idx = -1;
    std::unordered_set<int> seen_x;
    seen_x.reserve(pairs.size());
    nway_events.clear();
    nway_events.reserve(competitor_nodes.size());

    for (std::size_t i = 0; i < pairs.size(); ++i) {
      const Pair &p = pairs[i];
      bool l0_is_gate = (p.l0 == gate_label);
      bool l1_is_gate = (p.l1 == gate_label);
      if (l0_is_gate == l1_is_gate) {
        return false;
      }
      int gate_e = l0_is_gate ? p.e0 : p.e1;
      int x_e = l0_is_gate ? p.e1 : p.e0;
      int x_label = l0_is_gate ? p.l1 : p.l0;
      if (x_label == gate_label || x_label < 0) {
        return false;
      }
      if (!seen_x.insert(x_label).second) {
        return false;
      }
      if (i == 0) {
        gate_event_idx = gate_e;
        target_event_idx = x_e;
      } else {
        nway_events.push_back(x_e);
      }
    }
    if (gate_event_idx < 0 || target_event_idx < 0 || nway_events.empty()) {
      return false;
    }

    std::vector<int> event_node_indices;
    event_node_indices.reserve(2 + nway_events.size());
    event_node_indices.push_back(
        ir.events[static_cast<std::size_t>(gate_event_idx)].node_idx);
    event_node_indices.push_back(
        ir.events[static_cast<std::size_t>(target_event_idx)].node_idx);
    for (int ce : nway_events) {
      if (ce < 0 || ce >= static_cast<int>(ir.events.size())) {
        return false;
      }
      event_node_indices.push_back(
          ir.events[static_cast<std::size_t>(ce)].node_idx);
    }
    for (std::size_t i = 0; i < event_node_indices.size(); ++i) {
      int ni = event_node_indices[i];
      if (ni < 0 || ni >= static_cast<int>(ir.nodes.size())) {
        continue;
      }
      for (std::size_t j = i + 1; j < event_node_indices.size(); ++j) {
        int nj = event_node_indices[j];
        if (nj < 0 || nj >= static_cast<int>(ir.nodes.size())) {
          continue;
        }
        if (mask_overlap(ir.nodes[static_cast<std::size_t>(ni)],
                         ir.nodes[static_cast<std::size_t>(nj)])) {
          return false;
        }
      }
    }

    spec_out = IrSharedGateSpec();
    spec_out.kind = IrSharedGateKind::NWay;
    spec_out.node_idx = node_idx;
    spec_out.nway_gate_event_idx = gate_event_idx;
    spec_out.nway_target_event_idx = target_event_idx;
    spec_out.nway_competitor_event_count = static_cast<int>(nway_events.size());
    return true;
  };

  for (std::size_t oi = 0; oi < ir.outcomes.size(); ++oi) {
    IrOutcome &out = ir.outcomes[oi];
    if (out.node_idx < 0 || out.competitor_count <= 0) {
      continue;
    }
    std::vector<int> competitor_nodes;
    competitor_nodes.reserve(static_cast<std::size_t>(out.competitor_count));
    for (int k = 0; k < out.competitor_count; ++k) {
      competitor_nodes.push_back(ir.outcome_competitors[static_cast<std::size_t>(
          out.competitor_begin + k)]);
    }
    IrSharedGateSpec spec;
    bool ok = false;
    if (competitor_nodes.size() == 1) {
      ok = detect_shared_gate_pair_ir(out.node_idx, competitor_nodes[0], spec);
    } else if (competitor_nodes.size() >= 2) {
      std::vector<int> nway_events;
      ok = detect_shared_gate_nway_ir(out.node_idx, competitor_nodes, nway_events,
                                      spec);
      if (ok && !nway_events.empty()) {
        spec.nway_competitor_event_begin =
            static_cast<int>(ir.nway_competitor_event_indices.size());
        ir.nway_competitor_event_indices.insert(
            ir.nway_competitor_event_indices.end(), nway_events.begin(),
            nway_events.end());
      }
    }
    if (!ok) {
      continue;
    }
    spec.competitor_begin = out.competitor_begin;
    spec.competitor_count = out.competitor_count;
    int spec_idx = static_cast<int>(ir.shared_gate_specs.size());
    ir.shared_gate_specs.push_back(std::move(spec));
    out.shared_gate_spec_idx = spec_idx;
    std::vector<int> competitors;
    competitors.reserve(static_cast<std::size_t>(out.competitor_count));
    for (int k = 0; k < out.competitor_count; ++k) {
      competitors.push_back(ir.outcome_competitors[static_cast<std::size_t>(
          out.competitor_begin + k)]);
    }
    ir.shared_gate_lookup.emplace(shared_gate_lookup_key(out.node_idx, competitors),
                                  spec_idx);
  }

  ir.valid = true;
  ctx->ir = std::move(ir);
  ctx->kernel_program = compile_kernel_program(ctx->ir);
  ctx->base_params_soa = TrialParamsSoA{};
  ctx->base_params_soa.n_acc = static_cast<int>(ctx->accumulators.size());
  ctx->base_params_soa.dist_code.resize(ctx->accumulators.size(), 0);
  ctx->base_params_soa.onset.resize(ctx->accumulators.size(), 0.0);
  ctx->base_params_soa.q.resize(ctx->accumulators.size(), 0.0);
  ctx->base_params_soa.t0.resize(ctx->accumulators.size(), 0.0);
  ctx->base_params_soa.p1.resize(ctx->accumulators.size(), 0.0);
  ctx->base_params_soa.p2.resize(ctx->accumulators.size(), 0.0);
  ctx->base_params_soa.p3.resize(ctx->accumulators.size(), 0.0);
  for (std::size_t acc_i = 0; acc_i < ctx->accumulators.size(); ++acc_i) {
    const NativeAccumulator &acc = ctx->accumulators[acc_i];
    ctx->base_params_soa.dist_code[acc_i] = acc.dist_cfg.code;
    ctx->base_params_soa.onset[acc_i] = acc.onset;
    ctx->base_params_soa.q[acc_i] = acc.q;
    ctx->base_params_soa.t0[acc_i] = acc.dist_cfg.t0;
    ctx->base_params_soa.p1[acc_i] = acc.dist_cfg.p1;
    ctx->base_params_soa.p2[acc_i] = acc.dist_cfg.p2;
    ctx->base_params_soa.p3[acc_i] = acc.dist_cfg.p3;
  }
  ctx->base_params_soa.valid = true;
  ctx->kernel_state_graph = KernelStateGraph{};
  ctx->kernel_state_graph.forced_bit_count =
      static_cast<int>(ctx->ir.label_id_to_bit_idx.size());
  ctx->kernel_state_graph.bit_idx_to_label_id.assign(
      ctx->kernel_state_graph.forced_bit_count, NA_INTEGER);
  for (const auto &entry : ctx->ir.label_id_to_bit_idx) {
    const int label_id = entry.first;
    const int bit_idx = entry.second;
    if (bit_idx >= 0 &&
        bit_idx < static_cast<int>(ctx->kernel_state_graph.bit_idx_to_label_id.size())) {
      ctx->kernel_state_graph
          .bit_idx_to_label_id[static_cast<std::size_t>(bit_idx)] = label_id;
    }
  }
  ctx->kernel_state_graph.trigger_op_indices.clear();
  ctx->kernel_state_graph.trigger_transition_begin.assign(
      ctx->shared_trigger_groups.size(), -1);
  ctx->kernel_state_graph.trigger_transition_count.assign(
      ctx->shared_trigger_groups.size(), 0);
  ctx->kernel_state_graph.node_guard_transition_idx.assign(
      ctx->ir.nodes.size(), -1);
  ctx->kernel_state_graph.guard_transitions.clear();
  ctx->kernel_state_graph.guard_transitions.reserve(ctx->ir.nodes.size());
  for (int node_idx = 0; node_idx < static_cast<int>(ctx->ir.nodes.size());
       ++node_idx) {
    const IrNode &node = ctx->ir.nodes[static_cast<std::size_t>(node_idx)];
    if (node.source_mask_begin < 0 || node.source_mask_count <= 0) {
      continue;
    }
    KernelGuardTransition tr;
    tr.node_idx = node_idx;
    tr.source_mask_begin = node.source_mask_begin;
    tr.source_mask_count = node.source_mask_count;
    const int tr_idx =
        static_cast<int>(ctx->kernel_state_graph.guard_transitions.size());
    ctx->kernel_state_graph.guard_transitions.push_back(tr);
    ctx->kernel_state_graph.node_guard_transition_idx[static_cast<std::size_t>(
        node_idx)] = tr_idx;
  }
  std::size_t total_trigger_transitions = 0;
  for (const SharedTriggerGroup &group : ctx->shared_trigger_groups) {
    total_trigger_transitions += group.acc_indices.size();
  }
  ctx->kernel_state_graph.trigger_transitions.reserve(
      total_trigger_transitions);
  for (std::size_t tg = 0; tg < ctx->shared_trigger_groups.size(); ++tg) {
    const SharedTriggerGroup &group = ctx->shared_trigger_groups[tg];
    std::unordered_set<int> group_acc_indices(group.acc_indices.begin(),
                                              group.acc_indices.end());
    std::vector<int> affected_op_indices;
    affected_op_indices.reserve(ctx->kernel_program.ops.size());
    for (int op_idx = 0; op_idx < static_cast<int>(ctx->kernel_program.ops.size());
         ++op_idx) {
      const KernelOp &op = ctx->kernel_program.ops[static_cast<std::size_t>(op_idx)];
      if (op.code != KernelOpCode::Event) {
        continue;
      }
      if (op.event_idx < 0 ||
          op.event_idx >= static_cast<int>(ctx->ir.events.size())) {
        continue;
      }
      const IrEvent &event = ctx->ir.events[static_cast<std::size_t>(op.event_idx)];
      if (event.acc_idx < 0) {
        continue;
      }
      if (group_acc_indices.find(event.acc_idx) != group_acc_indices.end()) {
        affected_op_indices.push_back(op_idx);
      }
    }
    std::sort(affected_op_indices.begin(), affected_op_indices.end());
    affected_op_indices.erase(
        std::unique(affected_op_indices.begin(), affected_op_indices.end()),
        affected_op_indices.end());
    const int op_begin =
        static_cast<int>(ctx->kernel_state_graph.trigger_op_indices.size());
    ctx->kernel_state_graph.trigger_op_indices.insert(
        ctx->kernel_state_graph.trigger_op_indices.end(),
        affected_op_indices.begin(), affected_op_indices.end());
    const int op_count = static_cast<int>(affected_op_indices.size());
    const int transition_begin =
        static_cast<int>(ctx->kernel_state_graph.trigger_transitions.size());
    int transition_count = 0;
    for (int acc_idx : group.acc_indices) {
      KernelStateTransition tr;
      tr.trigger_bit = static_cast<int>(tg);
      tr.acc_idx = acc_idx;
      tr.op_begin = op_begin;
      tr.op_count = op_count;
      ctx->kernel_state_graph.trigger_transitions.push_back(tr);
      ++transition_count;
    }
    ctx->kernel_state_graph
        .trigger_transition_begin[static_cast<std::size_t>(tg)] =
        transition_begin;
    ctx->kernel_state_graph
        .trigger_transition_count[static_cast<std::size_t>(tg)] =
        transition_count;
  }
  validate_guard_transition_completeness(*ctx);
  ctx->kernel_state_graph.valid = ctx->kernel_program.valid;
  return ctx;
}

} // namespace uuber
