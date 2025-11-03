rm(list = ls())
source("R/generator_new.R"); 
source("examples/new_API.R"); 
source("R/super_large_likelihood.R"); 
examples <- names(new_api_examples)[c(1:11, 13:16, 18)]

for (nm in examples) {
  cat("Response Probabilities for model ", nm, ": \n")
  model <- new_api_examples[[nm]]; print(response_probabilities(model)) 
}
# Test output should be:
# Response Probabilities for model  example_1_simple  : 
#   R1        R2 
# 0.6000711 0.3999289 
# Response Probabilities for model  example_2_stop_mixture  : 
#   R1        R2 
# 0.8062219 0.1937782 
# Response Probabilities for model  example_3_stop_na  : 
#   Left     Right      STOP 
# 0.4068894 0.2657838 0.3273268 
# Response Probabilities for model  example_4_two_on_one  : 
#   R1        R2 
# 0.6640977 0.3359023 
# Response Probabilities for model  example_5_timeout_guess  : 
#   Left     Right   TIMEOUT 
# 0.4875719 0.5124281 0.0000000 
# Response Probabilities for model  example_6_dual_path  : 
#   Outcome_via_A Outcome_via_B 
# 0.7000557     0.2999443 
# Response Probabilities for model  example_7_mixture  : 
#   R1        R2 
# 0.3250216 0.6749784 
# Response Probabilities for model  example_8_shared_params  : 
#   Left Right 
# 0.5   0.5 
# Response Probabilities for model  example_9_advanced_k  : 
#   A         B 
# 0.2916013 0.7083987 
# Response Probabilities for model  example_10_exclusion  : 
#   R1        R2 
# 0.4678644 0.5321356 
# Response Probabilities for model  example_11_censor_deadline  : 
#   Left       Right   NR_CENSOR NR_DEADLINE 
# 0.7625428   0.2150445   0.0224127   0.0000000 
# Response Probabilities for model  example_12_inhibitor_with_protector  : 
#   R1        R2 
# 0.6970571 0.3029430 
# Response Probabilities for model  example_13_nested_pools  : 
#   TeamA     TeamB 
# 0.2956925 0.7043075 
# Response Probabilities for model  example_14_weighted_pool  : 
#   WeightedChoice     Competitor 
# 0.6568396      0.3431604 
# Response Probabilities for model  example_15_component_metadata  : 
#   Response     GUESS 
# 0.8199568 0.1800000 
# Response Probabilities for model  example_16_guard_tie_simple  : 
#   Fast      Slow 
# 0.4253897 0.5746103 
# Response Probabilities for model  example_18_shared_triggers  : 
#   Left     Right      STOP 
# 0.5123152 0.3470024 0.1406823 