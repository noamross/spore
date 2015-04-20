library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(noamtools)
library(nloptr)
library(polynom)
library(parallel)

closeAllConnections()
parms = list(
  max_i = 20,
  lambda = 0.001,
  lambda_ex = 0.2,
  alpha = 0.1,
  mu = 0.01,
  r = 0.5,
  d = 0.01,
  K = 100,
  init_pop = 100,
  time_max = 50,
  prevent_inf = 0,
  prevent_ex = 0,
  macro_timestep = 1,
  micro_timestep = 0.125,
  micro_relax_steps = 1,
  project = FALSE,
  n_sims = 100,
  control_min = 0,
  control_max = 1000,
  v = 10,
  c = 100,
  progress = FALSE
  #micro_record = file("micro.txt", open="w+")
  #  macro_record = file("macro.txt", open="w")
)


micro_state = c(100, rep(0, parms$max_i))
macro_state = restrict.micro_state(micro_state)
shadow_state = c(100, -100)
time = 0
# Rprof('opt.prof')
#a = determine_control(macro_state = macro_state, parms = parms, shadow_state = shadow_state, time = 0, control_guess = 0)
# Rprof(NULL)


options(mc.cores=2)
parms$control_max = 0

no_control_runs <- mclapply(1:50, function(x) macro_state_c_runopt(macro_state_init = macro_state, parms=parms, shadow_state_init=shadow_state, time=0, control_guess_init=0))

parms$control_max = 1000
control_runs <- mclapply(1:50, function(x) macro_state_c_runopt(macro_state_init = macro_state, parms=parms, shadow_state_init=shadow_state, time=0, control_guess_init=0))

parms$c = 10000
expensive_control_run = mclapply(1:50, function(x) (macro_state_c_runopt(macro_state_init = macro_state, parms=parms, shadow_state_init=shadow_state, time=0, control_guess_init=0))

