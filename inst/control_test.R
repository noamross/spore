devtools::load_all(".")
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(noamtools)
library(nloptr)
library(polynom)
library(parallel)
library(rlist)

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
  time_max = 2,
  prevent_inf = 0,
  prevent_ex = 0,
  macro_timestep = 0.25,
  micro_timestep = 0.025,
  micro_relax_steps = 3,
  project = FALSE,
  n_sims = 100,
  control_min = 0,
  control_max = 1000,
  v = 50,
  c = 100,
  progress = TRUE
  #micro_record = file("micro.txt", open="w+")
  #  macro_record = file("macro.txt", open="w")
)


micro_state = c(100, 0, rep(0, parms$max_i - 1))
macro_state = restrict.micro_state(micro_state)
shadow_state = c(100, -100)
time = 0
options(error=recover)
a = macro_state_c_runopt(macro_state_init = macro_state, parms=parms, shadow_state_init=shadow_state, time=0, control_guess_init=0)

#Rprof('opt.prof')
#for(i in 1:500) {
#a = determine_control(macro_state = macro_state, parms = parms, shadow_state = shadow_state, time = 0, control_guess = 0)
#}
# Rprof(NULL)


options(mc.cores=20)
#options(error = quote({dump.frames(to.file = TRUE)}))
parms$control_max = 0

no_control_runs <- mclapply(1:50, function(x) {
  myseed = sample.int(1e6,1)
  set.seed(myseed)
  out = try(macro_state_c_runopt(macro_state_init = macro_state, parms=parms, shadow_state_init=shadow_state, time=0, control_guess_init=0))
  if ("try-error" %in% class(out)) out = list(out, myseed)
  return(out)
})
saveRDS(no_control_runs, "no_control_runs.rds", compress=FALSE)

parms$control_max = 1000
control_runs <- mclapply(1:50, function(x) {
  myseed = sample.int(1e6, 1)
  set.seed(myseed)
  out = try(macro_state_c_runopt(macro_state_init = macro_state, parms=parms, shadow_state_init=shadow_state, time=0, control_guess_init=0))
  if ("try-error" %in% class(out)) out = list(out, myseed)
  return(out)
})
saveRDS(control_runs, "control_runs.rds", compress=FALSE)

parms$c = 10000
expensive_control_runs = mclapply(1:50, function(x) {
  myseed = sample.int(1e6,1)
  set.seed(myseed)
  out = try(macro_state_c_runopt(macro_state_init = macro_state, parms=parms, shadow_state_init=shadow_state, time=0, control_guess_init=0))
  if ("try-error" %in% class(out)) out = list(out, myseed)
  return(out)
})

saveRDS(expensive_control_runs, "expensive_control_runs.rds", compress=FALSE)


no_control_runs = readRDS('no_control_runs.rds')
control_runs = readRDS('control_runs.rds')
expensive_control_runs = readRDS('expensive_control_runs.rds')
process_runs = . %>%
  list.filter(!("try-error" %in% class(.))) %>%
  lapply(., as.data.frame) %>%
  list.map(cbind(run=.i, .)) %>%
  rbind_all %>%
  rename(run=run, time=times, N=V2, P=V3, dN=V4, dp=V5, ddN=V6, ddP=V7, S1 = V8, S2=V9, dS1=V10, dS2=V11, h = V12) %>%
  gather(key=variable, value=value, -run, - time)

no_control_runs_df = process_runs(no_control_runs)
control_runs_df = process_runs(control_runs)
expensive_control_runs_df = process_runs(expensive_control_runs)
library(ggplot2)

ggplot(subset(no_control_runs_df, variable %in% c("N", "P")), aes(x=time, y=value, col=variable, group=paste0(run,variable))) + geom_line(alpha = 0.5)
ggplot(subset(control_runs_df, variable %in% c("N", "P")), aes(x=time, y=value, col=variable, group=paste0(run,variable))) + geom_line(alpha = 0.5)
ggplot(subset(expensive_control_runs_df, variable %in% c("N", "P")), aes(x=time, y=value, col=variable, group=paste0(run,variable))) + geom_line(alpha = 0.5)
