rm(list = ls())
source("examples/stim_selective_versions.R")
source("R/model_tables.R")

nams <- names(new_api_examples)
for(nam in nams){
  cat(nam, "\n")
  print(model_to_tables(new_api_examples[[nam]]))
}

data.frame(trial = rep(1:4, each = 4), condition = rep(c("EASY", "HARD", "EASY", "HARD"), each = 4), 
           accumulator = rep(1:4, 4), type = rep(c("guard", "std"), 8),
           outcome = rep(c("R1", "R2", NA, "R1"), each = 4), rt = exp(rnorm(16)),
           meanlog = c(.1, .2, .3, .4), sdlog = c(.2, .2, .3, .3))