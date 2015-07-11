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
  max_i = 100,
  lambda = 0.001,
  lambda_ex = 0.2,
  alpha = 0.1,
  alpha_power = 1,
  mu = 0.01,
  r = 0.5,
  d = 0.01,
  K = 100,
  init_pop = 100,
  time_max = 40,
  prevent_inf = 0,
  prevent_ex = 0,
  macro_timestep = 0.5,
  micro_timestep = 0.05,
  micro_relax_steps = 1,
  project = FALSE,
  n_sims = 5000,
  control_min = 0,
  control_max = 1000,
  v = 50,
  c = 200,
  progress = TRUE
  #micro_record = file("micro.txt", open="w+")
  #  macro_record = file("macro.txt", open="w")
)


micro_state = c(100, 0, rep(0, parms$max_i - 1))
macro_state = restrict.micro_state(micro_state)
shadow_state = c(97.57118, -174.8631)
time = 0
#parms$control_max = 0

#Rprof("opt.prof")
b = macro_state_c_runopt(macro_state_init = macro_state, parms=parms, shadow_state_init=shadow_state, time=0, control_guess_init=0)
#Rprof(NULL)

process_runs = . %>%
  list.filter(!("try-error" %in% class(.))) %>%
  lapply(., as.data.frame) %>%
  list.map(cbind(run=.i, .)) %>%
  rbind_all %>%
  rename(run=run, time=times, N=V2, P=V3, S1=V4, S2=V5, dN=V6, dP=V7, dS1=V8,
         dS2=V9, ddNdN=V10, ddPdN=V11, ddNdP=V12, ddPdP=V13, h=V14,
         hamiltonian=hamiltonian) %>%
  gather(key=variable, value=value, -run, - time)

spread_runs = . %>% spread(variable, value) %>%
  group_by(run) %>%
  mutate(profit = sum(parms$v * N - parms$c*h)*parms$macro_timestep) %>%
  group_by()

run_profits = . %>% group_by(run) %>% summarize(value = sum(parms$v * N * parms$macro_timestep), cost = sum(parms$c *h * parms$macro_timestep), profit = value-cost)


a = process_runs(control_runs)
a2 = spread_runs(a)
profits = run_profits(a2)
ggplot(subset(a, variable %in% c("N", "P")), aes(x=time, y=value, col=variable, group=paste0(run,variable))) + geom_line(lwd=0.5)
ggplot(subset(a, variable %in% c("N", "P", "h")), aes(x=time, y=value, col=variable, group=paste0(run,variable))) + geom_line(lwd=0.5)
ggplot(subset(a, variable %in% c("S1", "S2") & time < 20), aes(x=time, y=value, col=variable, group=paste0(run,variable))) + geom_line(lwd=0.5)
ggplot(subset(a, variable %in% c("ddNdP")), aes(x=time, y=value, col=variable, group=paste0(run,variable))) + geom_line(lwd=0.5)


# sum(filter(a, variable=="N")$value*parms$v - filter(a, variable=="h")$value * parms$c)
profits %>% arrange(profit)

ggplot(subset(a2, time < 20), aes(x=time, y=h, group=run)) + geom_line()
hist(profits$profit)


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
saveRDS(no_control_runs, "no_control_runs_sh.rds", compress=FALSE)

parms$control_max = 1000
control_runs <- mclapply(1:20, function(x) {
  myseed = sample.int(1e6, 1)
  set.seed(myseed)
  out = try(macro_state_c_runopt(macro_state_init = macro_state, parms=parms, shadow_state_init=shadow_state, time=0, control_guess_init=0))
  if ("try-error" %in% class(out)) out = list(out, myseed)
  return(out)
})
saveRDS(control_runs, "control_runs_sh.rds", compress=FALSE)

parms$c = 10000
expensive_control_runs = mclapply(1:50, function(x) {
  myseed = sample.int(1e6,1)
  set.seed(myseed)
  out = try(macro_state_c_runopt(macro_state_init = macro_state, parms=parms, shadow_state_init=shadow_state, time=0, control_guess_init=0))
  if ("try-error" %in% class(out)) out = list(out, myseed)
  return(out)
})

saveRDS(expensive_control_runs, "expensive_control_runs_sh.rds", compress=FALSE)


ggplot(subset(a, run %in% 4 & variable %in% c("N", "P", "h")), aes(x=time, y=value, col=variable, group=paste0(run,variable))) + geom_line(lwd=0.5)
ggplot(subset(a, run %in% 4 & variable %in% c("S1", "S2") & time < 20), aes(x=time, y=value, col=variable, group=paste0(run,variable))) + geom_line(lwd=0.5)


ggplot(subset(a, run %in% 8:9 & variable %in% c("h")), aes(x=time, y=value, col=variable, group=paste0(run,variable))) + geom_line(lwd=0.5)
a2 %>% filter(run==5) %>% print(n=nrow(.))

