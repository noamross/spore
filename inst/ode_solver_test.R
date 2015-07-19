#devtools::load_all(".")
#library(spore)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(noamtools)
library(nloptr)
library(polynom)
library(parallel)
library(rlist)
library(tracer)
library(deSolve)

closeAllConnections()
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
  n_sims = 1000,
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


micro_state = c(1000, 0, rep(0, parms$max_i - 1))
macro_state = restrict.micro_state(micro_state)
names(macro_state) = c("N", "P")
shadow_state = c( 994.418947, -21566.367077 )
time = 0

ode_func = function(t, y, parms, shadow_state, control) {
  parms=relist(parms)
  a = macro_state_c_deriv_aves(macro_state=y, parms = parms, time = 0, control = control)
  output =  c(dN = a$macro_state_deriv[1],
               dP = a$macro_state_deriv[2])
  cat(c(y, output, t), "\n")
  return(list(output, output))
}
parmvec = unlist(as.relistable(parms))
out = ode(y = macro_state, times = 0:40, parms=parmvec, func = ode_func, method="adams", rtol = 1e-1, atol=1, shadow_state=shadow_state, control=0)

out_data = as.data.frame(out) %>%
  gather("variable", "value", -time)

init = c(macro_state, shadow_state); names(init) = c("N", "P", "S1", "S2")
file.remove('out.txt')
file.create('out.txt')
parmvec = unlist(as.relistable(parms))
system.time(out_opt <- ode(y = init, times = 0:40, parms=parmvec, func = opt_derivs, method = "euler"))
output = read.csv("out.txt", header=FALSE)
names(output) = c("time", "N", "P", "S1", "S2", "dN", "dP", "dS1", "dS2", "ddNdN", "ddPdN", "ddNdP", "ddPdP", "h", "H")

parms$n_sims = 100000

file.remove('out.txt')
file.create('out.txt')
parmvec = unlist(as.relistable(parms))
system.time(out_opt <- ode(y = init, times = 0:40, parms=parmvec, func = opt_derivs, method = "euler"))
output = read.csv("out.txt", header=FALSE)
names(output) = c("time", "N", "P", "S1", "S2", "dN", "dP", "dS1", "dS2", "ddNdN", "ddPdN", "ddNdP", "ddPdP", "h", "H")


