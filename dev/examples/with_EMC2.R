rm(list = ls())
library(EMC2)
library(AccumulatR)
load("~/Documents/2025/AccumulatR/dev/examples/test_AccumulatR.RData")

AccumulatR_add_ctx <- function(dadm, model_list){
  nacc <- unique
  data <- dadm[!duplicated(dadm$trial),]
  model_spec <- model_list$spec
  data$accumulator <- NULL
  ctx <- build_likelihood_context(model_spec, data)
  
  context <- list(
    native_ctx = ctx$native_ctx,  # externalptr (the “real” native context)
    rel_tol    = ctx$rel_tol,
    abs_tol    = ctx$abs_tol,
    max_depth  = ctx$max_depth
  )
  
  
  attr(dadm, "AccumulatR_context") <- context
  return(dadm)
}

names(finalize_model(model$spec)$prep$accumulators)

dadm <- AccumulatR_add_ctx(dadm, model)
context <- attr(dadm, "AccumulatR_context")
lls <- EMC2:::calc_ll_AccR(proposals, dadm, constants = constants, designs = designs,
                    model$bound, model$transform, model$pre_transform, p_types = p_types, min_ll = log(1e-10),
                    model$trend, context)
# rm(list = ls())
# library(EMC2)
# library(AccumulatR)
# load("~/Documents/2026/StopData/stop_opposite/dat_sc_emc.RData")
# dat_sc_emc$SS <- !is.infinite(dat_sc_emc$SSD)
# dat_sc_emc$SS <- factor(dat_sc_emc$SS)
# 
# plot_density(dat_sc_emc, factors = c("subjects", "S"), layout = c(1,2))
# 
# table(dat_sc_emc$R, dat_sc_emc$SS, dat_sc_emc$S)
# 
# dat <- dat_sc_emc[,c("subjects", "S", "R", "rt", "SSD", "SS")]
# dat[dat$subjects == unique(dat$subjects)[3],]
# dat$SSD[dat$SSD == Inf] <- 0
# dat <- dat[dat$rt > .15,]
# 
# dat <- dat[!is.na(dat$subjects),]
# dat$cond <- ifelse(dat$SSD == 0, "GO", "STOP")
# dat$cond <- factor(dat$cond)
# dat$component <- ifelse(dat$SSD == 0, "go_only", "go_stop")
# dat$component <- factor(dat$component)
# dat$subjects <- factor(dat$subjects)
# 
# 
# fail_go_idx <- dat$R == "Left" & dat$S == "Right" & dat$cond == "GO" |dat$R == "Right" & dat$S == "Left" & dat$cond == "GO"
# dat <- dat[!fail_go_idx,]
# 
# 
# univalent_stop_change <- race_spec() |>
#   add_accumulator("go_left", "exgauss") |>
#   add_accumulator("go_right", "exgauss") |>
#   add_accumulator("stop", "exgauss", onset = 0.15) |>
#   add_accumulator("change2left", "exgauss", onset = 0.20) |>
#   add_accumulator("change2right", "exgauss", onset = 0.20) |>
#   add_outcome("Left", first_of(inhibit("go_left", by = "stop"),all_of("stop", "change2left"))) |>
#   add_outcome("Right", first_of(inhibit("go_right", by = "stop"),all_of("stop", "change2right"))) |>
#   add_component("go_only", members = c("go_left", "go_right"), weight = .75) |>
#   add_component("go_stop",members = c("go_left", "go_right", "stop", "change2left", "change2right"), weight = .25) |>
#   set_mixture_options(mode = "fixed") |>
#   finalize_model()
# 
# stop_change_EMC2 <- EMC2:::AccumulatR_model(univalent_stop_change)
# 
# group_lmfun <- function(data){
#   out <- rep("go_match", nrow(data))
#   # First stop accumulator
#   idx_stop <- data$accumulator == "stop"
#   out[idx_stop] <- "stop"
#   # mismatching go accumulator
#   idx_go_mm <- data$accumulator == "go_left" & data$S == "Right" |
#     data$accumulator == "go_right" & data$S == "Left"
#   out[idx_go_mm] <- "go_mismatch"
#   # matching change accumulator
#   idx_change_m <- data$accumulator == "change2left" & data$S == "Right" |
#     data$accumulator == "change2right" & data$S == "Left"
#   out[idx_change_m] <- "change_match"
#   # mismatching change accumulator
#   idx_change_mm <- data$accumulator == "change2left" & data$S == "Left" |
#     data$accumulator == "change2right" & data$S == "Right"
#   out[idx_change_mm] <- "change_mismatch"
#   return(factor(out))
# }
# 
# lMfun <- function(data){
#   out <- rep(FALSE, nrow(data))
#   out[data$S == "Left" & (data$accumulator == "go_left" | data$accumulator == "change2right")] <- TRUE
#   out[data$S == "Right" & (data$accumulator == "go_right" | data$accumulator == "change2left")] <- TRUE
#   return(factor(out))
# }
# 
# group_fun <- function(data){
#   group <- rep("go", nrow(data))
#   idx_change <- data$accumulator %in% c("change2left", "change2right")
#   idx_stop <- data$accumulator == "stop"
#   group[idx_change] <- "change"
#   group[idx_stop] <- "stop"
#   return(factor(group))
# }
# 
# onset_fun <- function(data){
#   onset <- data$SSD
#   onset[data$accumulator %in% c("go_left", "go_right")] <- 0
#   return(onset)
# }
# 
# 
# ADmat <- cbind(d = c(1/2, -1/2))
# 
# des <- design(data = dat,
#               formula = list(p1 ~ 0 + gLM, p2 ~ 0 + gLM, p3 ~ 0 + group, t0 ~ 0 + group),
#               model = stop_change_EMC2,
#               functions = list(gLM = group_lmfun, group = group_fun, onset = onset_fun),
#               constants = c("t0_groupchange" = log(0),
#                             "t0_groupstop" = log(0),
#                             "p1_gLMgo_mismatch" = 100,
#                             "p1_gLMchange_mismatch" = 100,
#                             "p2_gLMgo_mismatch" = log(1),
#                             "p2_gLMchange_mismatch" = log(1)))
# emc <- make_emc(dat, des, type = "single")
# 
# debug(EMC2:::calc_ll_manager)
# emc <- fit(emc, fileName = "test_EMC2", iter = 250, cores_for_chains = 1)


# rm(list = ls())
# library(EMC2)
# library(AccumulatR)
# set.seed(123)
# R1 <- .1 + exp(rnorm(100, -.6, .45))
# R2 <- .1 + exp(rnorm(100, -.7, .3))
# 
# dat <- data.frame(rt = pmin(R1, R2), R = 1 + (R2 < R1))
# dat$R <- factor(dat$R)
# dat$subjects <- factor(1)
# 
# simple_exg <- race_spec() |>
#   add_accumulator("r1", "exgauss") |>
#   add_accumulator("r2", "exgauss") |>
#   add_outcome("1", "r1") |>
#   add_outcome("2", "r2") |>
#   finalize_model()
# 
# simple_exg_EMC2 <- EMC2:::AccumulatR_model(simple_exg)
# 
# 
# # debug(group_lmfun)
# des <- design(data = dat,
#               formula = list(p1 ~ accumulator, p2 ~ accumulator, p3 ~ accumulator, t0 ~ 1),
#               model = simple_exg_EMC2)
# 
# emc <- make_emc(dat, des, type = "single")
# emc <- fit(emc, fileName = "simple_exg", iter = 500)
# credint(emc)
# 
# # debug(make_data)
# pp <- predict(emc, n_cores = 10)
# 
# plot_cdf(dat, pp)
#
# credint(emc)
#
# data <- dat
# data$correct <- dat$S == dat$R
# table(data$correct, dat$cond)/nrow(dat)
#
# post <- pp
# post$correct <- pp$S == pp$R
# table(post$correct, post$cond)/nrow(post)
#
# credint(emc, map = TRUE)
# p_vector <- credint(emc, map = FALSE)[[1]][,2]
# dat <- make_data(p_vector, des, data = data)
# dat_rep <- rbind(dat, dat, dat, dat, dat, dat, dat)
# dat_rep$R <- factor(dat_rep$R)
# dat_rep$trial <- dat_rep$trials <- 1:nrow(dat_rep)
#
# devtools::load_all()
# profile_plot(dat, des, p_vector)
#
#
#
