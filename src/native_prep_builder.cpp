#include "native_prep_builder.hpp"

#include <algorithm>
#include <stdexcept>

#include "native_accumulator.hpp"

namespace uuber {
namespace {

std::vector<std::string> extract_string_vector(SEXP obj) {
  if (Rf_isNull(obj)) return {};
  Rcpp::CharacterVector chr(obj);
  std::vector<std::string> out;
  out.reserve(chr.size());
  for (R_xlen_t i = 0; i < chr.size(); ++i) {
    if (chr[i] == NA_STRING) continue;
    out.push_back(Rcpp::as<std::string>(chr[i]));
  }
  return out;
}

std::vector<int> extract_int_vector(SEXP obj, bool subtract_one = false) {
  if (Rf_isNull(obj)) return {};
  Rcpp::IntegerVector vec(obj);
  std::vector<int> out;
  out.reserve(vec.size());
  for (int val : vec) {
    if (val == NA_INTEGER) continue;
    out.push_back(subtract_one ? val - 1 : val);
  }
  return out;
}

std::string extract_string(SEXP obj) {
  if (Rf_isNull(obj) || obj == NA_STRING) return std::string();
  return Rcpp::as<std::string>(obj);
}

int extract_int(SEXP obj, int default_value = -1) {
  if (Rf_isNull(obj)) return default_value;
  int val = Rcpp::as<int>(obj);
  if (val == NA_INTEGER) return default_value;
  return val;
}

std::vector<ProtoParamEntry> extract_params(SEXP obj) {
  if (Rf_isNull(obj)) return {};
  Rcpp::List lst(obj);
  if (lst.size() == 0) return {};
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

NativePrepProto build_prep_proto(const Rcpp::List& prep) {
  NativePrepProto proto;

  Rcpp::List acc_list = prep["accumulators"];
  if (!acc_list.isNULL()) {
    Rcpp::CharacterVector acc_names(acc_list.names());
    for (R_xlen_t i = 0; i < acc_list.size(); ++i) {
      if (acc_list[i] == R_NilValue) continue;
      Rcpp::List acc(acc_list[i]);
      ProtoAccumulator acc_proto;
      if (!acc_names.isNULL() && acc_names[i] != NA_STRING) {
        acc_proto.id = Rcpp::as<std::string>(acc_names[i]);
      }
      if (acc.containsElementNamed("dist")) acc_proto.dist = extract_string(acc["dist"]);
      if (acc.containsElementNamed("onset") && !Rf_isNull(acc["onset"])) {
        acc_proto.onset = Rcpp::as<double>(acc["onset"]);
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
      if (pool_list[i] == R_NilValue) continue;
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
      if (names.isNULL() || names[i] == NA_STRING) continue;
      int value = extract_int(id_index[i], NA_INTEGER);
      if (value == NA_INTEGER) continue;
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
      if (node_list[i] == R_NilValue) continue;
      Rcpp::List node(node_list[i]);
      ProtoNode proto_node;
      proto_node.id = extract_int(node["id"], static_cast<int>(i + 1));
      if (node.containsElementNamed("kind")) proto_node.kind = extract_string(node["kind"]);
      if (node.containsElementNamed("source")) proto_node.source = extract_string(node["source"]);
      if (node.containsElementNamed("args")) proto_node.args = extract_int_vector(node["args"]);
      if (node.containsElementNamed("sources")) proto_node.source_ids = extract_int_vector(node["sources"]);
      if (node.containsElementNamed("reference_id")) proto_node.reference_id = extract_int(node["reference_id"]);
      if (node.containsElementNamed("blocker_id")) proto_node.blocker_id = extract_int(node["blocker_id"]);
      if (node.containsElementNamed("unless_ids")) proto_node.unless_ids = extract_int_vector(node["unless_ids"]);
      if (node.containsElementNamed("arg")) proto_node.arg_id = extract_int(node["arg"]);
      if (node.containsElementNamed("needs_forced")) {
        proto_node.needs_forced = Rcpp::as<bool>(node["needs_forced"]);
      }
      if (node.containsElementNamed("scenario_sensitive")) {
        proto_node.scenario_sensitive = Rcpp::as<bool>(node["scenario_sensitive"]);
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

  return proto;
}

std::unique_ptr<NativeContext> build_context_from_proto(const NativePrepProto& proto) {
  auto ctx = std::make_unique<NativeContext>();

  ctx->accumulators.reserve(proto.accumulators.size());
  for (const auto& acc_proto : proto.accumulators) {
    NativeAccumulator acc;
    acc.id = acc_proto.id;
    acc.dist = acc_proto.dist;
    acc.onset = acc_proto.onset;
    acc.q = acc_proto.q;
    acc.components = acc_proto.components;
    acc.shared_trigger_id = acc_proto.shared_trigger_id;
    try {
      acc.dist_cfg = resolve_acc_params_entries(acc.dist, acc_proto.params);
    } catch (const std::exception& e) {
      Rcpp::stop("native_context_build: accumulator '%s' %s",
                 acc.id.c_str(), e.what());
    }
    ctx->accumulators.push_back(acc);
  }

  for (std::size_t i = 0; i < ctx->accumulators.size(); ++i) {
    const auto& acc = ctx->accumulators[i];
    if (!acc.id.empty()) {
      ctx->accumulator_index[acc.id] = static_cast<int>(i);
    }
    if (!acc.shared_trigger_id.empty()) {
      ctx->shared_trigger_map[acc.shared_trigger_id].push_back(static_cast<int>(i));
    }
  }

  ctx->pools.reserve(proto.pools.size());
  for (const auto& pool_proto : proto.pools) {
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

  ctx->nodes.reserve(proto.nodes.size());
  for (const auto& node_proto : proto.nodes) {
    NativeNode node;
    node.id = node_proto.id;
    node.kind = node_proto.kind;
    node.source = node_proto.source;
    node.args = node_proto.args;
    node.source_ids = node_proto.source_ids;
    node.reference_id = node_proto.reference_id;
    node.blocker_id = node_proto.blocker_id;
    node.unless_ids = node_proto.unless_ids;
    node.arg_id = node_proto.arg_id;
    node.needs_forced = node_proto.needs_forced;
    node.scenario_sensitive = node_proto.scenario_sensitive;
    ctx->node_index[node.id] = static_cast<int>(ctx->nodes.size());
    ctx->nodes.push_back(std::move(node));
  }

  for (const auto& label : proto.label_index) {
    if (!label.label.empty()) {
      ctx->label_to_id[label.label] = label.id;
    }
  }

  return ctx;
}

} // namespace uuber
