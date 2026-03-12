#include "component_plan.h"

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

double normalize_weights(std::vector<double> &weights) {
  if (weights.empty()) {
    return 0.0;
  }
  bool need_uniform = false;
  double total = 0.0;
  for (double w : weights) {
    if (!std::isfinite(w) || w < 0.0) {
      need_uniform = true;
      break;
    }
    total += w;
  }
  if (!need_uniform && total > 0.0) {
    for (double &w : weights) {
      w /= total;
    }
    return total;
  }
  const double uniform = 1.0 / static_cast<double>(weights.size());
  std::fill(weights.begin(), weights.end(), uniform);
  return 1.0;
}

} // namespace

double extract_trial_id(const Rcpp::DataFrame *trial_rows) {
  if (!trial_rows) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (!trial_rows->containsElementNamed("trial")) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  Rcpp::NumericVector trial_col((*trial_rows)["trial"]);
  if (trial_col.size() == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const double val = trial_col[0];
  if (Rcpp::NumericVector::is_na(val)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return static_cast<double>(val);
}

bool build_component_mix(const Rcpp::CharacterVector &component_ids,
                         const Rcpp::NumericVector &weights_in,
                         Rcpp::Nullable<Rcpp::String> forced_component,
                         std::vector<std::string> &components_out,
                         std::vector<double> &weights_out) {
  if (component_ids.size() == 0) {
    return false;
  }
  components_out.clear();
  components_out.reserve(component_ids.size());
  for (R_xlen_t i = 0; i < component_ids.size(); ++i) {
    Rcpp::String comp(component_ids[i]);
    if (comp == NA_STRING) {
      components_out.emplace_back("__default__");
    } else {
      components_out.emplace_back(comp.get_cstring());
    }
  }
  weights_out.assign(components_out.size(),
                     std::numeric_limits<double>::quiet_NaN());
  if (weights_in.size() == component_ids.size()) {
    for (R_xlen_t i = 0; i < weights_in.size(); ++i) {
      weights_out[static_cast<std::size_t>(i)] = weights_in[i];
    }
  }
  normalize_weights(weights_out);
  if (forced_component.isNotNull()) {
    const std::string forced = Rcpp::as<std::string>(forced_component);
    auto it = std::find(components_out.begin(), components_out.end(), forced);
    if (it == components_out.end()) {
      return false;
    }
    const std::string selected = *it;
    components_out.clear();
    components_out.push_back(selected);
    weights_out.clear();
    weights_out.push_back(1.0);
  }
  return !components_out.empty();
}

Rcpp::List native_component_plan_impl(const Rcpp::List &structure,
                                      const Rcpp::DataFrame *trial_rows,
                                      double trial_id,
                                      const std::string *forced_component) {
  (void)trial_id;
  Rcpp::DataFrame comp_tbl(structure["components"]);
  Rcpp::CharacterVector component_ids = comp_tbl["component_id"];
  Rcpp::NumericVector base_weights = comp_tbl["weight"];
  std::string mode = "fixed";
  if (comp_tbl.containsElementNamed("mode")) {
    Rcpp::CharacterVector mode_vec(comp_tbl["mode"]);
    if (mode_vec.size() > 0 && mode_vec[0] != NA_STRING) {
      mode = Rcpp::as<std::string>(mode_vec[0]);
    }
  }
  std::string reference;
  if (comp_tbl.containsElementNamed("reference")) {
    Rcpp::CharacterVector ref_vec(comp_tbl["reference"]);
    if (ref_vec.size() > 0 && ref_vec[0] != NA_STRING) {
      reference = Rcpp::as<std::string>(ref_vec[0]);
    }
  }
  Rcpp::List attrs_list;
  if (comp_tbl.containsElementNamed("attrs")) {
    attrs_list = Rcpp::List(comp_tbl["attrs"]);
  }

  std::vector<std::string> all_components;
  all_components.reserve(component_ids.size());
  std::unordered_map<std::string, double> base_weight_map;
  base_weight_map.reserve(component_ids.size());
  for (R_xlen_t i = 0; i < component_ids.size(); ++i) {
    const std::string id = Rcpp::as<std::string>(component_ids[i]);
    const double weight =
        (i < base_weights.size()) ? static_cast<double>(base_weights[i]) : 0.0;
    all_components.push_back(id);
    base_weight_map[id] = weight;
  }

  std::vector<std::string> selected_components;
  if (forced_component && !forced_component->empty()) {
    selected_components.push_back(*forced_component);
  } else if (trial_rows) {
    Rcpp::CharacterVector trial_comp = (*trial_rows)["component"];
    std::unordered_set<std::string> trial_filter;
    trial_filter.reserve(static_cast<std::size_t>(trial_comp.size()));
    for (R_xlen_t i = 0; i < trial_comp.size(); ++i) {
      trial_filter.insert(Rcpp::as<std::string>(trial_comp[i]));
    }
    if (!trial_filter.empty()) {
      for (const auto &comp : all_components) {
        if (trial_filter.count(comp)) {
          selected_components.push_back(comp);
        }
      }
    }
  }
  if (selected_components.empty()) {
    selected_components = all_components;
  }

  std::vector<double> weights(selected_components.size(), 0.0);
  if (mode == "sample") {
    std::unordered_map<std::string, int> comp_index;
    for (std::size_t i = 0; i < selected_components.size(); ++i) {
      comp_index[selected_components[i]] = static_cast<int>(i);
    }
    const std::string ref_id =
        reference.empty() ? selected_components.front() : reference;
    double sum_nonref = 0.0;
    for (R_xlen_t i = 0; i < component_ids.size(); ++i) {
      const std::string cid = Rcpp::as<std::string>(component_ids[i]);
      if (cid == ref_id) {
        continue;
      }
      auto it_idx = comp_index.find(cid);
      if (it_idx == comp_index.end()) {
        continue;
      }
      std::string wp;
      if (attrs_list.size() > 0 && i < attrs_list.size()) {
        Rcpp::List attr(attrs_list[i]);
        if (attr.containsElementNamed("weight_param")) {
          wp = Rcpp::as<std::string>(attr["weight_param"]);
        }
      }
      if (wp.empty()) {
        Rcpp::stop("Component '%s' missing weight_param for sampled mixture",
                   cid.c_str());
      }
      double val = std::numeric_limits<double>::quiet_NaN();
      if (trial_rows && trial_rows->containsElementNamed(wp.c_str())) {
        Rcpp::NumericVector col((*trial_rows)[wp]);
        if (col.size() > 0) {
          const double cand = col[0];
          if (std::isfinite(cand)) {
            val = cand;
          }
        }
      }
      if (!std::isfinite(val)) {
        auto it = base_weight_map.find(cid);
        if (it != base_weight_map.end()) {
          val = it->second;
        }
      }
      if (!std::isfinite(val) || val < 0.0 || val > 1.0) {
        Rcpp::stop("Mixture weight '%s' must be a probability in [0,1]",
                   wp.c_str());
      }
      weights[static_cast<std::size_t>(it_idx->second)] = val;
      sum_nonref += val;
    }
    auto ref_it = comp_index.find(ref_id);
    if (ref_it != comp_index.end()) {
      double ref_weight = 1.0 - sum_nonref;
      if (ref_weight < -1e-8) {
        Rcpp::stop(
            "Mixture weights sum to >1; check non-reference weight params");
      }
      if (ref_weight < 0.0) {
        ref_weight = 0.0;
      }
      weights[static_cast<std::size_t>(ref_it->second)] = ref_weight;
    }
  } else {
    for (std::size_t i = 0; i < selected_components.size(); ++i) {
      auto it = base_weight_map.find(selected_components[i]);
      double w = (it != base_weight_map.end()) ? it->second : 0.0;
      if (!std::isfinite(w) || w < 0.0) {
        w = 0.0;
      }
      weights[i] = w;
    }
  }

  double total = 0.0;
  for (double w : weights) {
    total += w;
  }
  if (total <= 0.0) {
    const double uniform = 1.0 / static_cast<double>(weights.size());
    std::fill(weights.begin(), weights.end(), uniform);
  } else {
    for (double &w : weights) {
      w /= total;
    }
  }

  return Rcpp::List::create(
      Rcpp::Named("components") = Rcpp::CharacterVector(
          selected_components.begin(), selected_components.end()),
      Rcpp::Named("weights") =
          Rcpp::NumericVector(weights.begin(), weights.end()));
}
