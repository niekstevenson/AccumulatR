#pragma once

#include <Rcpp.h>

#include <string>
#include <vector>

double extract_trial_id(const Rcpp::DataFrame *trial_rows);
bool build_component_mix(const Rcpp::CharacterVector &component_ids,
                         const Rcpp::NumericVector &weights_in,
                         Rcpp::Nullable<Rcpp::String> forced_component,
                         std::vector<std::string> &components_out,
                         std::vector<double> &weights_out);
Rcpp::List native_component_plan_impl(const Rcpp::List &structure,
                                      const Rcpp::DataFrame *trial_rows,
                                      double trial_id,
                                      const std::string *forced_component);
