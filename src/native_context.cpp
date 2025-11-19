#include "native_context.h"

#include "native_prep_builder.hpp"

#include <limits>
#include <cmath>

namespace {

using namespace uuber;

void populate_component_metadata(const Rcpp::List& prep, NativeContext& ctx) {
  Rcpp::List components = prep["components"];
  if (components.isNULL()) return;
  Rcpp::CharacterVector comp_ids = components["ids"];
  Rcpp::List comp_attrs = components["attrs"];
  for (R_xlen_t i = 0; i < comp_ids.size(); ++i) {
    if (comp_ids[i] == NA_STRING) continue;
    std::string comp_label = Rcpp::as<std::string>(comp_ids[i]);
    ComponentContextInfo meta;
    meta.deadline = std::numeric_limits<double>::quiet_NaN();
    Rcpp::RObject attrs_obj = comp_attrs[i];
    if (!attrs_obj.isNULL()) {
      Rcpp::List attrs(attrs_obj);
      if (attrs.containsElementNamed("deadline")) {
        Rcpp::RObject deadline_obj = attrs["deadline"];
        if (!deadline_obj.isNULL()) {
          meta.deadline = Rcpp::as<double>(deadline_obj);
        }
      }
      if (attrs.containsElementNamed("guess")) {
        Rcpp::RObject guess_obj = attrs["guess"];
        if (!guess_obj.isNULL()) {
          Rcpp::List guess_list(guess_obj);
          if (guess_list.containsElementNamed("outcome")) {
            meta.guess.target = Rcpp::as<std::string>(guess_list["outcome"]);
          }
          if (meta.guess.target.empty()) {
            meta.guess.target = "__GUESS__";
          }
          if (guess_list.containsElementNamed("weights")) {
            Rcpp::NumericVector weights(guess_list["weights"]);
            Rcpp::CharacterVector weight_names = weights.names();
            if (!weight_names.isNULL()) {
              for (R_xlen_t j = 0; j < weights.size(); ++j) {
                if (weight_names[j] == NA_STRING) continue;
                std::string source_label = Rcpp::as<std::string>(weight_names[j]);
                double keep_prob = weights[j];
                meta.guess.keep_weights[source_label] = keep_prob;
              }
            }
          }
          meta.guess.valid = !meta.guess.keep_weights.empty();
        }
      }
    }
    ctx.component_info[comp_label] = meta;
  }
}

void populate_outcome_metadata(const Rcpp::List& prep, NativeContext& ctx) {
  Rcpp::List outcomes = prep["outcomes"];
  if (!outcomes.isNULL()) {
    Rcpp::CharacterVector outcome_names = outcomes.names();
    for (R_xlen_t i = 0; i < outcomes.size(); ++i) {
      if (outcome_names[i] == NA_STRING) continue;
      Rcpp::RObject def_obj = outcomes[i];
      if (def_obj.isNULL()) continue;
      std::string outcome_label = Rcpp::as<std::string>(outcome_names[i]);
      OutcomeContextInfo info;
      info.node_id = -1;
      Rcpp::List def(def_obj);
      if (def.containsElementNamed("expr")) {
        Rcpp::RObject expr = def["expr"];
        if (!expr.isNULL()) {
          Rcpp::RObject lik_id = expr.attr(".lik_id");
          if (!lik_id.isNULL()) {
            info.node_id = Rcpp::as<int>(lik_id);
          }
        }
      }
      Rcpp::List options = def.containsElementNamed("options") ? Rcpp::List(def["options"]) : Rcpp::List();
      if (!options.isNULL()) {
        if (options.containsElementNamed("map_outcome_to")) {
          Rcpp::RObject map_obj = options["map_outcome_to"];
          if (!map_obj.isNULL()) {
            if (map_obj == NA_STRING) {
              ctx.alias_sources["__NA__"].push_back(outcome_label);
              info.maps_to_na = true;
            } else if (Rf_isString(map_obj)) {
              std::string target = Rcpp::as<std::string>(map_obj);
              ctx.alias_sources[target].push_back(outcome_label);
            }
          }
        }
        if (options.containsElementNamed("guess")) {
          Rcpp::RObject guess_obj = options["guess"];
          if (!guess_obj.isNULL()) {
            Rcpp::List guess_list(guess_obj);
            Rcpp::CharacterVector labels = guess_list.containsElementNamed("labels")
              ? Rcpp::CharacterVector(guess_list["labels"])
              : Rcpp::CharacterVector();
            Rcpp::NumericVector weights = guess_list.containsElementNamed("weights")
              ? Rcpp::NumericVector(guess_list["weights"])
              : Rcpp::NumericVector();
            std::string rt_policy = "keep";
            if (guess_list.containsElementNamed("rt_policy")) {
              Rcpp::RObject pol = guess_list["rt_policy"];
              if (!pol.isNULL()) {
                rt_policy = Rcpp::as<std::string>(pol);
              }
            }
            if (labels.size() == weights.size() && labels.size() > 0) {
              for (R_xlen_t j = 0; j < labels.size(); ++j) {
                if (labels[j] == NA_STRING) continue;
                OutcomeGuessDonor donor;
                donor.label = outcome_label;
                donor.weight = weights[j];
                donor.rt_policy = rt_policy;
                std::string target = Rcpp::as<std::string>(labels[j]);
                ctx.guess_target_map[target].push_back(std::move(donor));
              }
            }
          }
        }
      }
      ctx.outcome_info[outcome_label] = info;
    }
  }

  Rcpp::List competitors = prep[".competitors"];
  if (!competitors.isNULL()) {
    Rcpp::CharacterVector comp_names = competitors.names();
    for (R_xlen_t i = 0; i < competitors.size(); ++i) {
      if (comp_names[i] == NA_STRING) continue;
      Rcpp::RObject comp_obj = competitors[i];
      if (comp_obj.isNULL()) continue;
      std::string outcome_label = Rcpp::as<std::string>(comp_names[i]);
      Rcpp::List comp_list(comp_obj);
      std::vector<int> comp_ids;
      comp_ids.reserve(comp_list.size());
      for (R_xlen_t j = 0; j < comp_list.size(); ++j) {
        Rcpp::RObject entry = comp_list[j];
        if (entry.isNULL()) continue;
        Rcpp::RObject lik_id = entry.attr(".lik_id");
        if (lik_id.isNULL()) continue;
        comp_ids.push_back(Rcpp::as<int>(lik_id));
      }
      auto it = ctx.outcome_info.find(outcome_label);
      if (it == ctx.outcome_info.end()) {
        OutcomeContextInfo info;
        info.competitor_ids = comp_ids;
        ctx.outcome_info[outcome_label] = info;
      } else {
        it->second.competitor_ids = comp_ids;
      }
    }
  }
}

int resolve_na_cache_limit() {
  int default_limit = 128;
  Rcpp::Function getOption("getOption");
  Rcpp::RObject opt = getOption("uuber.cache.na.max_per_trial", Rcpp::wrap(default_limit));
  if (opt.isNULL()) return default_limit;
  try {
    double val = Rcpp::as<double>(opt);
    if (!std::isfinite(val) || val < 0.0) return 0;
    return static_cast<int>(std::floor(val));
  } catch (...) {
    return default_limit;
  }
}

} // namespace

namespace uuber {

Rcpp::XPtr<NativeContext> build_native_context(Rcpp::List prep) {
  NativePrepProto proto = build_prep_proto(prep);
  std::unique_ptr<NativeContext> ctx = build_context_from_proto(proto);
  populate_component_metadata(prep, *ctx);
  populate_outcome_metadata(prep, *ctx);
  ctx->na_cache_limit = resolve_na_cache_limit();
  ctx->context_builds = 1;
  ctx->context_reuses = 0;
  return Rcpp::XPtr<NativeContext>(ctx.release(), true);
}

} // namespace uuber
