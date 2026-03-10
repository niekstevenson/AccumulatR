#pragma once

#include <vector>

double pool_density_fast(const std::vector<double> &density,
                         const std::vector<double> &survival, int k);

double pool_survival_fast(const std::vector<double> &survival, int k);
