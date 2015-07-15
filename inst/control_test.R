#devtools::load_all(".")
library(spore)
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
  project = FALSE,
  n_sims = 100,
  control_min = 0,
  control_max = 0,
  v = 500,
  c = 2000,
  progress = TRUE,
  micro_record = 0,
  parallel_cores = 23
  #micro_record = file("micro.txt", open="w+")
  #  macro_record = file("macro.txt", open="w")
)


micro_state = c(1000, 0, rep(0, parms$max_i - 1))
macro_state = restrict.micro_state(micro_state)
shadow_state = c( 994.418947, -21566.367077 )
time = 0
#parms$control_max = 0
#determine_control(macro_state = macro_state, parms = parms, shadow_state = shadow_state, time = time, control_guess = 3, verbose = TRUE)
#Rprof("opt.prof")
b = macro_state_c_runopt(macro_state_init = macro_state, parms=parms, shadow_state_init=shadow_state, time=0, control_guess_init=3)
#Rprof(NULL)

noamtools::proftable('opt.prof')

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


a = process_runs(list(b))
a2 = spread_runs(a)
profits = run_profits(a2)
ggplot(subset(a, variable %in% c("N", "P")), aes(x=time, y=value, col=variable, group=paste0(run,variable))) + geom_line(lwd=0.5)
ggplot(subset(a, variable %in% c("N", "P", "h")), aes(x=time, y=value, col=variable, group=paste0(run,variable))) + geom_line(lwd=0.5)
ggplot(subset(a, variable %in% c("S1", "S2")), aes(x=time, y=value, col=variable, group=paste0(run,variable))) + geom_line(lwd=0.5)
ggplot(subset(a, variable %in% c("ddNdN", "ddPdN", "ddNdP", "ddPdP")), aes(x=time, y=value, col=variable, group=paste0(run,variable))) + geom_line(lwd=0.5)

run_means = a %>%
  group_by(variable, time) %>%
  summarize(value = mean(value))

ggplot(subset(run_means, variable %in% c("ddNdN", "ddPdN", "ddNdP", "ddPdP")), aes(x=time, y=value, col=variable)) + geom_line(lwd=0.5)


# sum(filter(a, variable=="N")$value*parms$v - filter(a, variable=="h")$value * parms$c)
profits %>% arrange(profit)

ggplot(subset(a2, time < 20), aes(x=time, y=h, group=run)) + geom_line()
hist(profits$profit)


options(mc.cores=20)
#options(error = quote({dump.frames(to.file = TRUE)}))
parms$control_max = 0

no_control_runs <- mclapply(1:20, function(x) {
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

parms$parallel_cores = 1
#a = plyr::rlply(100, sims_step(macro_state, parms, shadow_state, 0, 2.7), .progress="time", .parallel = TRUE)
parms$n_sims = 500
a1 = mclapply(1:1000, function(i) {sims_step(macro_state, parms, shadow_state, 0, 2.7)}, mc.preschedule = TRUE,  mc.cores = 23)
parms$n_sims = 2500
a2 = mclapply(1:1000, function(i) {sims_step(macro_state, parms, shadow_state, 0, 2.7)}, mc.preschedule = TRUE, mc.cores = 23)
parms$n_sims = 10000
a3 = mclapply(1:1000, function(i) {sims_step(macro_state, parms, shadow_state, 0, 2.7)}, mc.preschedule = TRUE, mc.cores = 23)

run_data = . %>%
  list.map(as.list(c(.$macro_state_deriv, as.vector(.$macro_dfdx)))) %>%
  list.stack %>%
  rename(dN = V1, dP = V2, ddNdN = V3, ddPdN = V4, ddNdP = V5, ddPdP = V6) %>%
  mutate(run = 1:n()) %>%
  gather("variable", "value", -run) %>%
  mutate(vartype = ifelse(variable %in% c("dN", "dP"), "derivatives", "partials"))

a1r = run_data(a1); a1r$nsims = 500;
a2r = run_data(a2); a2r$nsims = 2500;
a3r = run_data(a3); a3r$nsims = 10000;

a_data = rbind(a1r, a2r, a3r)

ggplot(a_data, aes(x=value, fill = variable, col=variable)) + geom_density(alpha = 0.5) + facet_grid(nsims~vartype, scales = "free")

parms$parallel_cores = 23
c = determine_control(macro_state, parms, shadow_state, time, 2.7, TRUE, TRUE, nloptr_options = list(algorithm = "NLOPT_LN_BOBYQA", xtol_rel = 1e-3))
d = determine_control(macro_state, parms, shadow_state, time, 2.7, FALSE, TRUE, nloptr_options = list(algorithm = "NLOPT_LN_SBPLX", xtol_rel = 1e-3))
c$opt$iterations
d$opt$iterations
plot(tracer(c$opt)$xval1, tracer(c$opt)$fval)
plot(tracer(d$opt)$xval1, tracer(d$opt)$fval)


plot(tracer(c$opt)$xval1, tracer(c$opt)$fval, xlim=c(2.5, 3.2), ylim=c(-4200, -4000))
plot(tracer(d$opt)$xval1, tracer(d$opt)$fval, xlim=c(2.5, 3.2), ylim=c(-4200, -4000))

system.time(sims_step(macro_state, parms, shadow_state, 0, 2.7))


parms = list(
  max_i = 100,
  is = 0:100,
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
  macro_timestep = 1,
  micro_timestep = 0.01,
  micro_relax_steps = 0,
  project = FALSE,
  n_sims = 10000,
  control_min = 0,
  control_max = 100,
  v = 50,
  c = 200,
  progress = TRUE,
  micro_record = 0,
  parallel_cores = 23
  #micro_record = file("micro.txt", open="w+")
  #  macro_record = file("macro.txt", open="w")
)

parms$n_sims = 500
e1 = mclapply(1:1000, function(i) {sims_step(macro_state, parms, shadow_state, 0, 2.7)}, mc.cores = 24)
parms$n_sims = 2500
e2 = mclapply(1:1000, function(i) {sims_step(macro_state, parms, shadow_state, 0, 2.7)}, mc.cores = 24)
parms$n_sims = 10000
e3 = mclapply(1:1000, function(i) {sims_step(macro_state, parms, shadow_state, 0, 2.7)}, mc.cores = 24)

run_data = . %>%
  list.map(as.list(c(.$macro_state_deriv, as.vector(.$macro_dfdx)))) %>%
  list.stack %>%
  rename(dN = V1, dP = V2, ddNdN = V3, ddPdN = V4, ddNdP = V5, ddPdP = V6) %>%
  mutate(run = 1:n()) %>%
  gather("variable", "value", -run) %>%
  mutate(vartype = ifelse(variable %in% c("dN", "dP"), "derivatives", "partials"))

e1r = run_data(a1); a1r$nsims = 500;
e2r = run_data(a2); a2r$nsims = 2500;
e3r = run_data(a3); a3r$nsims = 10000;

e_data = rbind(e1r, e2r, e3r)

ggplot(e_data, aes(x=value, fill = variable, col=variable)) + geom_density(alpha = 0.5) + facet_grid(nsims~vartype, scales = "free")

