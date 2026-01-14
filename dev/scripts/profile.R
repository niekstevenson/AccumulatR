rm(list = ls())
library(EMC2)
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
  add_accumulator("go_left", "exgauss") |>
  add_accumulator("go_right", "exgauss") |>
  add_accumulator("stop", "exgauss") |>
  add_accumulator("change", "exgauss") |>
  add_pool("GO_LEFT", "go_left") |>
  add_pool("GO_RIGHT", "go_right") |>
  add_pool("STOP", "stop") |>
  add_pool("CHANGE", "change") |>
  add_outcome("A", inhibit("GO_LEFT", by = "STOP"),
              options = list(component = c("go_only", "go_stop"))) |>
  add_outcome("B", inhibit("GO_RIGHT", by = "STOP"),
              options = list(component = c("go_only", "go_stop"))) |>
  add_outcome("X", all_of("CHANGE", "STOP"),
              options = list(component = "go_stop")) |>
  add_group("component:go_only", members = c("go_left", "go_right"),
            attrs = list(component = "go_only")) |>
  add_group("component:go_stop",
            members = c("go_left", "go_right", "stop", "change"),
            attrs = list(component = "go_stop")) |>
  set_metadata(mixture = list(mode = "fixed", # component weights are only relevant for data simulation
                              components = list(component("go_only", weight = 0.75),
                                                component("go_stop", weight = 0.25))
  )) |>
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
              formula = list(p1 ~ 1, p2 ~ 1, t0 ~ 1, p3 ~ 1),
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
proposals <- rmvnorm(n_particles, mean = c(-.5, log(.5), log(.2), log(.2)), sigma = diag(.1, 4))
colnames(proposals) <- names(EMC2::sampled_pars(emc))

system.time(
  lls <- EMC2:::calc_ll_AccR(proposals, dadm, constants = constants, designs = designs,
                             model$bound, model$transform, model$pre_transform, p_types = p_types, min_ll = log(1e-10),
                             model$trend, context)
)




