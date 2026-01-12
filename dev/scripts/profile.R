rm(list = ls())
library(EMC2)
load("~/Documents/2026/StopData/profile.RData")
load("~/Documents/2026/StopData/stop_change_3B.RData")
emc <- EMC2:::AccumulatR_check_context(emc)
dadm <- emc[[1]]$data[[1]]
constants <- attr(dadm, "constants")
context <- attr(dadm, "AccumulatR_context")

lls <- EMC2:::calc_ll_AccR(proposals, dadm, constants = constants, designs = designs,
                             model$bound, model$transform, model$pre_transform, p_types = p_types, min_ll = log(1e-10),
                             model$trend, context)


