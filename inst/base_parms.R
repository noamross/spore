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
library(nloptr)
library(R.cache)
library(extrafont)

## ----setparms------------------------------------------------------------
parms = list(
  max_i = 100,
  is = 0:100,
  lambda = 0.0001,
  lambda_ex = 0.5,
  alpha = 0.2,
  alpha_power = 1,
  mu = 0.01,
  r = 0.2,
  d = 0.01,
  K = 500,
  init_pop = 300,
  time_max = 10,
  prevent_inf = 0,
  prevent_ex = 0,
  macro_timestep = 1,
  micro_timestep = 0.1,
  micro_relax_steps = 0,
  delta = 0,
  project = FALSE,
  n_comp_sims = 100,
  n_sims = 1000,
  n_sims_jacob = 10000000,
  control_min = 0,
  control_max = 10,
  v = 4,
  c = 400,
  nocontrol = TRUE,
  progress = TRUE,
  micro_record = 0,
  parallel_cores = 2
  #micro_record = file("micro.txt", open="w+")
  #  macro_record = file("macro.txt", open="w")
)

micro_state = c(295, 5, rep(0, parms$max_i - 1))

macro_state = restrict.micro_state(micro_state)

N0 = 299.9999
P0 = 5.00001
init = c(N = N0, P = P0, S1 = 0, S2 = 0)
times = seq(0,parms$time_max, by=1)

