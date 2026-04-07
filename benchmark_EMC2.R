rm(list = ls())
library(EMC2)
library(AccumulatR)
load("~/Documents/2025/EMC2/test_accumulatR.RData")

likelihood_context <- make_context(model$spec)
dadm <- prepare_data(model$spec, as.data.frame(dadm))
p_types <- names(model$p_types)
lls <- EMC2:::calc_ll_oo(proposals, dadm, constants = constants, designs = designs, type = model$c_name,
                  model$bound, model$transform, model$pre_transform, p_types = p_types,
                  min_ll = log(1e-10), trend = model$trend,
                  likelihood_context = likelihood_context)
# rm(list = ls())
# # devtools::load_all()
# d <- read.csv("~/Documents/2026/StopData/3B_stopchange/stop_change_v1_pilot_2_data.csv")
# library(EMC2)
# library(AccumulatR)
#
#
#
# # Some inspection
# table(d$R)
# table(d$R)/nrow(d)
# table(d$S)
# head(d)
#
#
# # Some d wrangling
# d <- d[,c("s", "SS", "S", "R", "RT", "SSD")]
# colnames(d)[c(1,5)] <- c("subjects", "rt")
# d$SSD[d$SS == "GoTask"] <- 0
#
# d <- d[d$subjects %in% c("8953", "8954"),]
#
#
# fail_go_idx <- d$R == "A" & d$S == "B" & d$SS == "GoTask" |d$R == "B" & d$S == "A" & d$SS == "GoTask"
# mean(fail_go_idx)
# d <- d[!fail_go_idx,]
#
# # Wrong x
# fail_x_idx <- d$R == "X" & d$SS == "GoTask"
# mean(fail_x_idx)
# d <- d[!fail_x_idx,]
#
# d$S <- "S"
# d$R <- as.character(d$R)
# d$R[d$R != "X"] <- "S"
#
# range(d$rt)
#
# table(d$SS, d$S, d$R)
#
# # Recode levels
# d$S <- factor(d$S)
# d$SS <- factor(d$SS)
# d$R <- factor(d$R)
# d$subjects <- factor(d$subjects)
#
# d$component <- "go_only"
# d$component[d$SS == "ChangeTask"] <- "go_stop"
# d$component <- factor(d$component)
#
# model <- race_spec() |>
#   add_accumulator("S", "lognormal") |>
#   add_accumulator("stop", "lognormal") |>
#   add_accumulator("change", "lognormal") |>
#   add_outcome("S", inhibit("S", by = "stop")) |>
#   add_outcome("X", all_of("change", "stop")) |>
#   add_component("go_only", members = "S", weight = .75) |>
#   add_component("go_stop", members = c("S", "stop", "change"), weight = .25) |>
#   add_trigger(
#     "stop_trigger",
#     members = c("stop", "change"),
#     q = 0.05,
#     param = "stop_trigger"
#   ) |>
#   set_mixture_options(mode = "fixed") |>
#   finalize_model()
#
#
#
# stop_change_EMC2 <- EMC2:::AccumulatR_model(model)
#
#
# onset_fun <- function(data){
#   onset <- data$SSD
#   onset[data$accumulator %in% "S"] <- 0
#   return(onset)
# }
#
# type_fun <- function(data){
#   out <- rep("GO", nrow(data))
#   out[data$accumulator %in% c("stop", "change")] <- "STOP"
#   return(factor(out))
# }
#
#
#
# fit_subject <- function(i, data){
#   dat <- data[data$subjects==unique(data$subjects)[i],]
#   dat$subjects <- factor(dat$subjects)
#   des <- design(data = dat,
#                 formula = list(p1 ~ 0 + accumulator, p2 ~ 0 + accumulator, t0 ~ 0 + accumulator,
#                                q ~ 0 + type),
#                 model = stop_change_EMC2,
#                 functions = list(onset = onset_fun, type = type_fun),
#                 constants = c("q_typeGO" = qnorm(0),
#                               "t0_accumulatorstop" = log(0))
#   )
#   filename <- paste0(i, "3B_Spair_lnr_fixed.RData")
#   emc <- make_emc(dat, des, type = "single", compress = T)
#   emc <- fit(emc, fileName = filename, iter = 250, cores_for_chains = 1)
#   return(emc)
# }
#
# library(parallel)
# debug(EMC2:::calc_ll_manager)
# emcs <- mclapply(1:length(unique(d$subjects)), fit_subject, d, mc.cores = 1)
