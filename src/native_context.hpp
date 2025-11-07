#pragma once

#include "native_context.h"

namespace uuber {
namespace {

std::vector<std::string> extract_string_vector(SEXP obj) {
  if (Rf_isNull(obj)) return {};
  Rcpp::CharacterVector chr;
  try {
    chr = Rcpp::as<Rcpp::CharacterVector>(obj);
  } catch (const std::exception& e) {
    Rcpp::stop("native_context_build: expected character vector, got type %s (%s)",
               Rf_type2char(TYPEOF(obj)), e.what());
  }
  std::vector<std::string> out;
  out.reserve(chr.size());
  for (R_xlen_t i = 0; i < chr.size(); ++i) {
    if (chr[i] == NA_STRING) continue;
    out.push_back(Rcpp::as<std::string>(chr[i]));
  }
  return out;
}

std::vector<int> extract_int_vector(SEXP obj, bool subtract_one = false) {
  std::vector<int> out;
  if (Rf_isNull(obj)) return out;
  Rcpp::IntegerVector vec(obj);
  out.reserve(vec.size());
  for (int val : vec) {
    if (val == NA_INTEGER) continue;
    if (subtract_one) {
      out.push_back(val - 1);
    } else {
      out.push_back(val);
    }
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

template <typename T>
std::unordered_map<std::string, int> build_index(const std::vector<T>& items) {
  std::unordered_map<std::string, int> index;
  for (std::size_t i = 0; i < items.size(); ++i) {
    if (!items[i].id.empty()) {
      index[items[i].id] = static_cast<int>(i);
    }
  }
  return index;
}

} // namespace

inline Rcpp::XPtr<NativeContext> build_native_context(Rcpp::List prep) {
  auto ctx = std::make_unique<NativeContext>();

  Rcpp::List acc_list = prep["accumulators"];
  if (!acc_list.isNULL()) {
    Rcpp::CharacterVector acc_names;
    SEXP acc_names_sexp = acc_list.names();
    if (!Rf_isNull(acc_names_sexp)) {
      acc_names = Rcpp::CharacterVector(acc_names_sexp);
    }
    for (R_xlen_t i = 0; i < acc_list.size(); ++i) {
      try {
        Rcpp::List acc = acc_list[i];
        NativeAccumulator acc_def;
        if (!acc_names.isNULL() && acc_names[i] != NA_STRING) {
          acc_def.id = extract_string(acc_names[i]);
        }
        if (acc.containsElementNamed("dist")) {
          SEXP dist = acc["dist"];
          acc_def.dist = extract_string(dist);
        }
        if (acc.containsElementNamed("onset")) {
          SEXP onset = acc["onset"];
          if (!Rf_isNull(onset)) {
            acc_def.onset = Rcpp::as<double>(onset);
          }
        }
        if (acc.containsElementNamed("q")) {
          SEXP q = acc["q"];
          if (!Rf_isNull(q)) {
            acc_def.q = Rcpp::as<double>(q);
          }
        }
        if (acc.containsElementNamed("params")) {
          SEXP params = acc["params"];
          if (!Rf_isNull(params)) {
            acc_def.params = Rcpp::List(params);
          }
        } else {
          acc_def.params = Rcpp::List();
        }
        if (acc.containsElementNamed("components")) {
          SEXP comps = acc["components"];
          if (!Rf_isNull(comps)) {
            acc_def.components = extract_string_vector(comps);
          }
        }
        if (acc.containsElementNamed("shared_trigger_id")) {
          SEXP sid = acc["shared_trigger_id"];
          acc_def.shared_trigger_id = extract_string(sid);
        }
        ctx->accumulators.push_back(std::move(acc_def));
      } catch (const std::exception& e) {
        Rcpp::stop("native_context_build: accumulator[%d] %s", static_cast<int>(i + 1), e.what());
      }
    }
  }
  Rcpp::List pool_list = prep["pools"];
  if (!pool_list.isNULL()) {
    Rcpp::CharacterVector pool_names;
    SEXP pool_names_sexp = pool_list.names();
    if (!Rf_isNull(pool_names_sexp)) {
      pool_names = Rcpp::CharacterVector(pool_names_sexp);
    }
    for (R_xlen_t i = 0; i < pool_list.size(); ++i) {
      try {
        Rcpp::List pool = pool_list[i];
        NativePool pool_def;
        if (!pool_names.isNULL() && pool_names[i] != NA_STRING) {
          pool_def.id = extract_string(pool_names[i]);
        }
        if (pool.containsElementNamed("k")) {
          SEXP kval = pool["k"];
          if (!Rf_isNull(kval)) {
            pool_def.k = Rcpp::as<int>(kval);
          }
        }
        if (pool.containsElementNamed("members")) {
          SEXP members = pool["members"];
          if (!Rf_isNull(members)) {
            pool_def.members = extract_string_vector(members);
          }
        }
        ctx->pools.push_back(std::move(pool_def));
      } catch (const std::exception& e) {
        Rcpp::stop("native_context_build: pool[%d] %s", static_cast<int>(i + 1), e.what());
      }
    }
  }
  ctx->accumulator_index = build_index(ctx->accumulators);
  ctx->pool_index = build_index(ctx->pools);

  Rcpp::List id_index = prep[".id_index"];
  if (!id_index.isNULL()) {
    Rcpp::CharacterVector names;
    SEXP names_sexp = id_index.names();
    if (!Rf_isNull(names_sexp)) {
      names = Rcpp::CharacterVector(names_sexp);
    }
    for (R_xlen_t i = 0; i < id_index.size(); ++i) {
      if (names.isNULL() || names[i] == NA_STRING) continue;
      int value = extract_int(id_index[i], NA_INTEGER);
      if (value == NA_INTEGER) continue;
      ctx->label_to_id[extract_string(names[i])] = value;
    }
  }

  Rcpp::List expr_compiled = prep[".expr_compiled"];
  if (!expr_compiled.isNULL()) {
    Rcpp::List node_list = expr_compiled["nodes"];
    ctx->nodes.reserve(node_list.size());
    for (R_xlen_t i = 0; i < node_list.size(); ++i) {
      if (node_list[i] == R_NilValue) continue;
      try {
        Rcpp::List node = node_list[i];
        NativeNode native;
        native.id = extract_int(node["id"], static_cast<int>(i + 1));
        if (node.containsElementNamed("kind")) {
          native.kind = extract_string(node["kind"]);
        }
        if (node.containsElementNamed("args")) {
          native.args = extract_int_vector(node["args"]);
        }
        if (node.containsElementNamed("sources")) {
          native.source_ids = extract_int_vector(node["sources"]);
        }
        if (node.containsElementNamed("reference_id")) {
          native.reference_id = extract_int(node["reference_id"]);
        }
        if (node.containsElementNamed("blocker_id")) {
          native.blocker_id = extract_int(node["blocker_id"]);
        }
        if (node.containsElementNamed("unless_ids")) {
          native.unless_ids = extract_int_vector(node["unless_ids"]);
        }
        if (node.containsElementNamed("arg")) {
          native.arg_id = extract_int(node["arg"]);
        }
        if (node.containsElementNamed("needs_forced")) {
          native.needs_forced = Rcpp::as<bool>(node["needs_forced"]);
        }
        if (node.containsElementNamed("scenario_sensitive")) {
          native.scenario_sensitive = Rcpp::as<bool>(node["scenario_sensitive"]);
        }
        if (node.containsElementNamed("expr")) {
          Rcpp::List expr = node["expr"];
          if (!expr.isNULL() && expr.containsElementNamed("source")) {
            native.source = extract_string(expr["source"]);
          }
        }
        ctx->node_index[native.id] = static_cast<int>(ctx->nodes.size());
        ctx->nodes.push_back(std::move(native));
      } catch (const std::exception& e) {
        Rcpp::stop("native_context_build: node[%d] %s", static_cast<int>(i + 1), e.what());
      }
    }
  }
  return Rcpp::XPtr<NativeContext>(ctx.release(), true);
}

} // namespace uuber
