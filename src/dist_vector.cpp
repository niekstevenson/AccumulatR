#include "dist_vector.h"

#include <algorithm>
#include <cstddef>

#include "accumulator.h"
#include "distribution_kernels.h"
#include "native_utils.h"

namespace uuber {
namespace {

template <typename Kernel>
inline void fill_pdf_vector(const double *x, std::size_t n, double *out,
                            Kernel kernel) {
  for (std::size_t i = 0; i < n; ++i) {
    out[i] = safe_density(kernel(x[i]));
  }
}

template <typename Kernel>
inline void fill_cdf_vector(const double *x, std::size_t n, double *out,
                            Kernel kernel) {
  for (std::size_t i = 0; i < n; ++i) {
    out[i] = clamp_probability(kernel(x[i]));
  }
}

} // namespace

void eval_pdf_vec(int dist_code, double p1, double p2, double p3, double p4,
                  double p5, double p6, double p7, double p8,
                  const double *x, std::size_t n, double *out) {
  if (!x || !out) {
    return;
  }
  (void)p5;
  (void)p6;
  (void)p7;
  (void)p8;
  switch (dist_code) {
  case ACC_DIST_LOGNORMAL:
    fill_pdf_vector(x, n, out, [=](double xi) {
      return uuber::distkernels::lognormal_pdf(xi, p1, p2);
    });
    return;
  case ACC_DIST_GAMMA:
    fill_pdf_vector(x, n, out, [=](double xi) {
      return uuber::distkernels::gamma_pdf(xi, p1, p2);
    });
    return;
  case ACC_DIST_EXGAUSS:
    fill_pdf_vector(x, n, out, [=](double xi) {
      return uuber::distkernels::exgauss_pdf(xi, p1, p2, p3);
    });
    return;
  case ACC_DIST_LBA:
    fill_pdf_vector(x, n, out, [=](double xi) {
      return uuber::distkernels::lba_pdf(xi, p1, p2, p3, p4);
    });
    return;
  case ACC_DIST_RDM:
    fill_pdf_vector(x, n, out, [=](double xi) {
      return uuber::distkernels::rdm_pdf(xi, p1, p2, p3, p4);
    });
    return;
  default:
    std::fill(out, out + n, 0.0);
    return;
  }
}

void eval_cdf_vec(int dist_code, double p1, double p2, double p3, double p4,
                  double p5, double p6, double p7, double p8,
                  const double *x, std::size_t n, double *out) {
  if (!x || !out) {
    return;
  }
  (void)p5;
  (void)p6;
  (void)p7;
  (void)p8;
  switch (dist_code) {
  case ACC_DIST_LOGNORMAL:
    fill_cdf_vector(x, n, out, [=](double xi) {
      return uuber::distkernels::lognormal_cdf(xi, p1, p2);
    });
    return;
  case ACC_DIST_GAMMA:
    fill_cdf_vector(x, n, out, [=](double xi) {
      return uuber::distkernels::gamma_cdf(xi, p1, p2);
    });
    return;
  case ACC_DIST_EXGAUSS:
    fill_cdf_vector(x, n, out, [=](double xi) {
      return uuber::distkernels::exgauss_cdf(xi, p1, p2, p3);
    });
    return;
  case ACC_DIST_LBA:
    fill_cdf_vector(x, n, out, [=](double xi) {
      return uuber::distkernels::lba_cdf(xi, p1, p2, p3, p4);
    });
    return;
  case ACC_DIST_RDM:
    fill_cdf_vector(x, n, out, [=](double xi) {
      return uuber::distkernels::rdm_cdf(xi, p1, p2, p3, p4);
    });
    return;
  default:
    std::fill(out, out + n, 0.0);
    return;
  }
}

} // namespace uuber
