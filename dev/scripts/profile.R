rm(list = ls())
devtools::load_all()
library(AccumulatR)
set.seed(123)

# Supplies data
load("~/Documents/2026/StopData/3B_stopchange/clean_data_stop_change_v1_pilot.Rdata")

data <- d

n_particles <- 100

data <- data[,c("s", "SS", "S", "R", "RT", "SSD")]
colnames(data)[c(1,5)] <- c("subjects", "rt")
data <- data[data$subjects == unique(data$subjects)[1],]

# Recode levels
data$S <- factor(data$S)
data$R <- factor(data$R)
data$subjects <- factor(data$subjects)
data$SSD[data$SS == "GoTask"] <- 0


model <- race_spec() |>
  add_accumulator("go_left", "lognormal") |>
  add_accumulator("go_right", "lognormal") |>
  add_accumulator("stop", "lognormal") |>
  add_accumulator("change", "lognormal") |>
  add_outcome("A", inhibit("go_left", by = "stop"),
              options = list(component = c("go_only", "go_stop"))) |>
  add_outcome("B", inhibit("go_right", by = "stop"),
              options = list(component = c("go_only", "go_stop"))) |>
  add_outcome("X", all_of("change", "stop"),
              options = list(component = "go_stop")) |>
  add_component("go_only", members = c("go_left", "go_right"), weight = .75) |>
  add_component("go_stop", members = c("go_left", "go_right", "stop", "change"), weight = .25) |>
  set_mixture_options(mode = "fixed") |>
  add_trigger(
    "shared_trigger",
    members = c("stop", "change"),
    q = 0.10,
    param = "stop_q",
    draw = "shared"
  ) |>
  finalize_model()


to_EMC2 <- EMC2:::AccumulatR_model(model)

onset_fun <- function(data){
  onset <- data$SSD
  onset[data$accumulator %in% c("go_left", "go_right")] <- 0
  return(onset)
}

data$component <- "go_only"
data$component[data$SS == "ChangeTask"] <- "go_stop"
data$component <- factor(data$component)

des <- design(data = data,
              formula = list(p1 ~ 1, p2 ~ 1, t0 ~ 1, q ~ 1),
              functions = list(onset = onset_fun),
              model = to_EMC2)

emc <- make_emc(data, des, type = "single")
model <- emc[[1]]$model()
dadm <- emc[[1]]$data[[1]]

p_types <- names(model$p_types)
designs <- list()
for(p in p_types){
  designs[[p]] <- attr(dadm,"designs")[[p]][attr(attr(dadm,"designs")[[p]],"expand"),,drop=FALSE]
}
constants <- attr(dadm, "constants")
context <- attr(dadm, "AccumulatR_context")

library(mvtnorm)
proposals <- rmvnorm(n_particles, mean = c(-.5, log(.5), log(.2), qnorm(.1)), sigma = diag(.1, 4))
colnames(proposals) <- names(EMC2::sampled_pars(emc))

system.time(
  lls <- EMC2:::calc_ll_AccR(proposals, dadm, constants = constants, designs = designs,
                             model$bound, model$transform, model$pre_transform, p_types = p_types, min_ll = log(1e-10),
                             model$trend, context)
)




