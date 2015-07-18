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
  n_sims = 1,
  control_min = 0,
  control_max = 0,
  v = 500,
  c = 2000,
  nocontrol = FALSE,
  progress = TRUE,
  micro_record = 0,
  parallel_cores = 3
  #micro_record = file("micro.txt", open="w+")
  #  macro_record = file("macro.txt", open="w")
)


N0 = 1000
P0 = 0
S10 =  100
S20 =  -400
init = c(N = N0, P = P0, S1 = S10, S2 = S20)
times = seq(0,40, by=1)

time = times[1]
micro_state = c(1000, 0, rep(0, parms$max_i - 1))
macro_state = restrict.micro_state(micro_state)
shadow_state = c( 994.418947, -21566.367077 )

constctrl = function(time) return(0)


library(parallel)
simulate_macro_state = function(macro_state, parms, times, controlfn, n_runs) {
  runs = mclapply(1:n_runs, mc.preschedule = TRUE, mc.cores = parms$parallel_cores,
                  FUN=function(run) {
                    output = data.frame(run = factor(run), time = times, N = NA, P = NA)
                    output[1, 3:4] = macro_state
                    macro_state_old = macro_state
                    for(i in 2:length(times)) {
                      macro_state_new = macro_state_c_stepto(macro_state_old, parms, controlfn(times[i]), times[i-1], times[i])
                      output[i, 3:4] = macro_state_new
                      macro_state_old = macro_state_new
                    }
                    return(output)
                  })
  return(dplyr::bind_rows(runs))
}

library(dplyr)
library(tidyr)
sim_runs = simulate_macro_state(macro_state, parms, times, constctrl, 1000) %>%
  gather(variable, value, N, P)
sim_run_subset = sim_runs %>% filter(run %in% sample.int(max(as.integer(run)), 10))
sim_aves = sim_runs %>%
  group_by(time, variable) %>%
  summarize(mean = mean(value), lower_sim=quantile(value, 0.05), upper_sim = quantile(value, 0.95),
            std_err = sd(value)/sqrt(length(value)), lower = mean - std_err, upper = mean + std_err)

h_fun = function(x) {
  ifelse(x >= 1, log(x), 0)
}

syseq2 = function(t, state, parms) {
  with(as.list(c(state, parms)), {
    #    h_calc = c(parms$control_min, parms$control_max, uniroot.all(hH_h,c(parms$control_min, parms$control_max), n=10000, t=t, state=state, parms=parms))
    #    if(length(h_calc) == 3) h_calc = c(h_calc, NA)
    #    H_calc = H(h_calc, t, state, parms)
    #    if(length(H_calc) == 3) H_calc = c(h_calc, NA)
    #    H_sel = max(H_calc, na.rm=TRUE)
    #    h = h_calc[which.max(H_calc)]
    if(parms$nocontrol) {
      h = 0
    } else {
      h = h_fun(-S2 * N * lambda_ex / c)
    }
    #h = suppressWarnings()
    #if(is.nan(h)) h = 0
    #h = min(h, parms$control_max)
    dN = max(r*N*(1 - N/K), 0) - alpha*P - d*N
    if(isTRUE(all.equal(dN, 0))) dN = 0
    dP = lambda*P*N - mu*P - d*P - alpha*P - alpha*(P^2)/N + exp(-h)*N*lambda_ex
    ddNdN= r - d - 2*(r/K)*N
    ddPdN=lambda*P + alpha*(P^2/N^2) + exp(-h) * lambda_ex
    ddNdP=-alpha
    ddPdP=lambda*N - mu - d - alpha - 2*alpha*P/N
    dS1 = -v - S1*ddNdN - S2*ddPdN + delta*S1
    dS2 = -S1*ddNdP - S2*ddPdP + delta*S2
    derivs = c(dN=dN, dP=dP, dS1=dS1, dS2=dS2)
    return(list(derivs, c(derivs, h=h, ddNdN=ddNdN, ddPdN=ddPdN, ddNdP=ddNdP, ddPdP=ddPdP)))
  })
}


parms$nocontrol = TRUE
library(deSolve)
init = c(N = macro_state[1], P = macro_state[2], S1=0, S2=0)
ode_no_control = as.data.frame(lsoda(init, times, syseq2, parms))
ode_no_control_data = ode_no_control %>%
  select(time, N, P) %>%
  gather(variable, value, -time)

init = c(N = macro_state[1], P = macro_state[2], S1=994.418947, S2=-21566.367077)
parms$nocontrol = FALSE
ode_control = as.data.frame(lsoda(init, times, syseq2, parms))
ode_control_data = ode_control %>%
  select(time, N, P, h) %>%
  gather(variable, value, -time)

controlfun = approxfun(ode_control$h)

library(ggplot2)
library(noamtools)

sim_runs_control = simulate_macro_state(macro_state, parms, times, controlfun, 1000) %>%
  group_by(run) %>%
  arrange(time) %>%
  mutate(h = ode_control$h) %>%
  group_by() %>%
  gather(variable, value, N, P)

sim_run_control_subset = sim_runs_control %>% filter(run %in% sample.int(max(as.integer(run)), 10))
sim_aves_control = sim_runs_control %>%
  group_by(time, variable) %>%
  summarize(mean = mean(value), lower_sim=quantile(value, 0.05), upper_sim = quantile(value, 0.95),
            std_err = sd(value)/sqrt(length(value)), lower = mean - std_err, upper = mean + std_err)

ggplot() +
  geom_ribbon(data = sim_aves_control, mapping = aes(x = time, ymin = lower_sim, ymax = upper_sim, group=variable),
              fill = "grey", color = "darkgrey", alpha = 0.5) +
  geom_line(data = sim_run_control_subset, mapping = aes(x = time, y = value, color = variable, group = paste(run, variable)), alpha = 1, lwd=0.25) +
  geom_ribbon(data = sim_aves_control, mapping = aes(x = time, ymin = lower, ymax = upper, fill=variable),
               alpha = 1) +
  geom_line(data = sim_aves_control, mapping = aes(x = time, y = mean, group = variable), color = "black", alpha = 1, lwd=0.5) +
 # geom_line(data = ode_control_data, mapping = aes(x = time, y = value, color = variable), alpha = 1, lty = 2, lwd = 1.5) +
  theme_nr

profout = function(out) {
  intval = integrate(approxfun(out$time, (out$N*parms$v - (out$h)*parms$c)*exp(-out$time*parms$delta)),
            lower = out$time[1], upper=tail(out$time, 1))
  return(intval$value)
}

ode_profit = profout(ode_control)
profits = sim_runs_control %>%
  spread(variable, value) %>%
  group_by(run) %>%
  do({
    data.frame(profit = profout(.))
  })

ggplot() +
  geom_density(data = profits, aes(x=profit), fill="green", col= NA, alpha = 0.4) +
  geom_vline(xintercept = ode_profit)
  xlim(19e6, 19.4e6) +
  theme_nr

parms$n_sims = 1000
eq = macro_state_c_runopt(macro_state_init = macro_state, parms=parms, shadow_state_init=shadow_state, time=0, control_guess_init=0)
library(rlist)
process_runs = . %>%
  list.filter(!("try-error" %in% class(.))) %>%
  lapply(., as.data.frame) %>%
  list.map(cbind(run=.i, .)) %>%
  rbind_all %>%
  rename(run=run, time=times, N=V2, P=V3, S1=V4, S2=V5, dN=V6, dP=V7, dS1=V8,
         dS2=V9, ddNdN=V10, ddPdN=V11, ddNdP=V12, ddPdP=V13, h=V14,
         hamiltonian=hamiltonian) %>%
  gather(key=variable, value=value, -run, - time)
eq_data = process_runs(list(eq))


ggplot() +
#   geom_ribbon(data = sim_aves, mapping = aes(x = time, ymin = lower_sim, ymax = upper_sim, group=variable),
#               fill = "grey", color = "darkgrey", alpha = 0.5) +
#   geom_line(data = sim_run_subset, mapping = aes(x = time, y = value, color = variable, group = paste(run, variable)), alpha = 1, lwd=0.25) +
#   geom_ribbon(data = sim_aves, mapping = aes(x = time, ymin = lower, ymax = upper, fill=variable),
#                alpha = 1) +
   geom_line(data = sim_aves, mapping = aes(x = time, y = mean, group = variable), color = "black", alpha = 1, lwd=0.5) +
   geom_line(data = ode_no_control_data, mapping = aes(x = time, y = value, color = variable), alpha = 1, lty = 1, lwd = 0.5) +
  geom_line(data = subset(eq_data, variable %in% c("N", "P")), mapping = aes(x = time, y = value, color = variable), lty=1, lwd = 0.5) +
  theme_nr

  ggplot() +
    geom_line(data = ode_no_control_data, mapping = aes(x = time, y = value, color = variable), alpha = 1, lty = 2, lwd = 1.5) +
    geom_line(data = ode_control_data, mapping = aes(x = time, y = value, color = variable), alpha = 1, lty = 1, lwd = 1.5) +
    theme_nr
