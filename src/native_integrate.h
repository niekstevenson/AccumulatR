#pragma once

#include <Rcpp.h>

namespace uuber {

double integrate_boost(Rcpp::Function integrand,
                       double lower,
                       double upper,
                       double rel_tol,
                       double abs_tol,
                       int max_depth);

} // namespace uuber
