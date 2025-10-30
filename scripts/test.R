rm(list = ls())
source("R/generator_new.R"); 
source("examples/new_API.R"); 
source("R/super_large_likelihood.R"); 
examples <- c("example_2_stop_mixture","example_6_dual_path","example_16_guard_tie_simple"); 
for (nm in examples) {
  model <- new_api_examples[[nm]]; print(response_probabilities(model)) 
}
# Outcome should be:
# R1        R2 
# 0.8062219 0.1937782 
# Outcome_via_A Outcome_via_B 
# 0.7000557     0.2999443 
# Fast      Slow 
# 0.4253897 0.5746103 