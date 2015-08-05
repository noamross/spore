#This chunk sets the base parameters for all of the simulations.

## ----loadpkgs------------------------------------------------------------

#library(spore)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(scales)
library(zoo)
library(spore)
library(readr)
library(deSolve)
library(parallel)
library(noamtools)
library(cowplot)


## ----setparms------------------------------------------------------------
parms = list(
  max_i = 100,
  is = 0:100,
  lambda = 0.0002,
  lambda_ex = 0.05,
  alpha = 0.1,
  alpha_power = 1,
  mu = 0.01,
  r = 0.5,
  d = 0.01,
  K = 1000,
  init_pop = 1000,
  time_max = 40,
  prevent_inf = 0,
  prevent_ex = 0,
  macro_timestep = 1,
  micro_timestep = 0.1,
  micro_relax_steps = 0,
  delta = 0,
  project = FALSE,
  n_comp_sims = 100,
  n_sims = 10000,
  n_sims_jacob = 10000000,
  control_min = 0,
  control_max = 10,
  v = 500,
  c = 2000,
  nocontrol = TRUE,
  progress = TRUE,
  micro_record = 0,
  parallel_cores = 23
  #micro_record = file("micro.txt", open="w+")
  #  macro_record = file("macro.txt", open="w")
)

micro_state = c(1000, rep(0, parms$max_i))

macro_state = restrict.micro_state(micro_state)

