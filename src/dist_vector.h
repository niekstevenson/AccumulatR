#pragma once

#include <cstddef>

namespace uuber {

void eval_pdf_vec(int dist_code, double p1, double p2, double p3, double p4,
                  double p5, double p6, double p7, double p8,
                  const double *x, std::size_t n, double *out);

void eval_cdf_vec(int dist_code, double p1, double p2, double p3, double p4,
                  double p5, double p6, double p7, double p8,
                  const double *x, std::size_t n, double *out);

} // namespace uuber
