// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>

#include <cstdlib>
#include <string>

#include <gperftools/profiler.h>

namespace {

void set_profile_frequency(const int frequency) {
  if (frequency <= 0) {
    return;
  }
  const auto value = std::to_string(frequency);
  setenv("CPUPROFILE_FREQUENCY", value.c_str(), 1);
}

} // namespace

// [[Rcpp::export]]
bool semantic_profiler_is_available_cpp() {
  return true;
}

// [[Rcpp::export]]
void semantic_profiler_start_cpp(SEXP pathSEXP, SEXP frequencySEXP) {
  const auto path = Rcpp::as<std::string>(pathSEXP);
  const auto frequency = Rcpp::as<int>(frequencySEXP);
  if (path.empty()) {
    Rcpp::stop("profile path must not be empty");
  }
  if (ProfilingIsEnabledForAllThreads() != 0) {
    ProfilerStop();
  }
  set_profile_frequency(frequency);
  if (ProfilerStart(path.c_str()) == 0) {
    Rcpp::stop("ProfilerStart failed");
  }
}

// [[Rcpp::export]]
void semantic_profiler_flush_cpp() {
  ProfilerFlush();
}

// [[Rcpp::export]]
void semantic_profiler_stop_cpp() {
  if (ProfilingIsEnabledForAllThreads() != 0) {
    ProfilerStop();
  }
}
