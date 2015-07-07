library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)

closeAllConnections()
parms = list(
  max_i = 20,
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
  micro_timestep = 0.1,
  micro_relax_steps = 3,
  project = FALSE,
  n_sims = 2000,
  control_min = 0,
  control_max = 1000,
  v = 50,
  c = 200,
  progress = TRUE
  #micro_record = file("micro.txt", open="w+")
  #  macro_record = file("macro.txt", open="w")
)

micro_state = c(100, rep(0, parms$max_i))

macro_state = restrict.micro_state(micro_state)

for(run in 1:50) {
	micro_state.stepto(micro_state, parms = parms, time = 0, timeto=parms$time_max, record = parms$micro_record, run = run)
}
close(parms$micro_record)

micro_output = load_micro_output('micro.txt')

Rprof('restrict_data.prof')
macro_output_short = micro_output %>%
  arrange(run, start, time) %>%
	group_by(run, start, time) %>%
  summarise(N = sum(population), P = sum(infections*population)) %>%
	mutate(time = plyr::round_any(time, 1, f=floor)) %>%
	filter(!duplicated(time)) %>%
	group_by()
Rprof(NULL)

macro_output_mean = macro_output_short %>%
	group_by(time) %>%
	summarize(mean_N = mean(N), mean_P = mean(P))

ggplot() +
	geom_line(data = macro_output_short, mapping = aes(x = time, y = N, group = run), col="green", alpha = 0.1, lwd=1) +
	geom_line(data = macro_output_short, mapping = aes(x = time, y = P, group = run), col="red", alpha = 0.1, lwd=1) +
	geom_line(data = macro_output_mean, mapping = aes(x = time, y = mean_N), col="green", lwd = 1) +
	geom_line(data = macro_output_mean, mapping = aes(x = time, y = mean_P), col="red", lwd = 1)


parms$micro_record = file("micro2.txt", open="w+")


macro_out = macro_state.run(macro_state, parms, time = 0)
close(parms$micro_record)

macro_out = tbl_df(as.data.frame(macro_out))
names(macro_out) = c("time", "N", "P")


micro_output_stepfull = load_micro_output('micro2.txt')
macro_output_stepfull = micro_output_stepfull %>%
  arrange(start, time, run) %>%
  group_by(run, start, time) %>%
  summarise(N = sum(population), P = sum(infections*population))


ggplot() +
	geom_line(data = macro_output_stepfull, mapping = aes(x=time, y = N, group = paste(run, start)), col = "green", alpha = 0.5) +
#	geom_line(data = macro_output_stepfull, mapping = aes(x=time, y = P, group = paste(run, start)), col = "red", alpha = 1) +
	geom_line(data = macro_out, mapping = aes(x = time, y = N), col="green", lwd=1) +
	geom_point(data = macro_out, mapping = aes(x = time, y = N)) +
#  geom_line(data = macro_out, mapping = aes(x = time, y = P), col="red", lwd=1) +
	xlim(0, 0.5) + ylim(90, 100)

ggplot() +
	geom_line(data = macro_output_short, mapping = aes(x = time, y = N, group = run), col="green", alpha = 0.1, lwd=1) +
	geom_line(data = macro_output_short, mapping = aes(x = time, y = P, group = run), col="red", alpha = 0.1, lwd=1) +
	geom_line(data = macro_output_mean, mapping = aes(x = time, y = mean_N), col="green", lwd = 1) +
	geom_line(data = macro_output_mean, mapping = aes(x = time, y = mean_P), col="red", lwd = 1) +
	geom_line(data = macro_out, mapping = aes(x = time, y = N), col="green", lwd=1, linetype = 2) +
	geom_line(data = macro_out, mapping = aes(x = time, y = P), col="red", lwd=1, linetype = 2)

parms$micro_record = file("micro3.txt", open="w+")
parms$project = TRUE
macro_proj_out = macro_state.run(macro_state, parms, time = 0)
close(parms$micro_record)

macro_proj_out = tbl_df(as.data.frame(macro_proj_out))
names(macro_proj_out) = c("time", "N", "P")

micro_output_stepproj = load_micro_output('micro3.txt')
macro_output_stepproj = micro_output_stepproj %>%
  arrange(start, time, run) %>%
  group_by(run, start, time) %>%
  summarise(N = sum(population), P = sum(infections*population))


ggplot() +
	geom_line(data = macro_output_short, mapping = aes(x = time, y = N, group = run), col="green", alpha = 0.1, lwd=1) +
	geom_line(data = macro_output_short, mapping = aes(x = time, y = P, group = run), col="red", alpha = 0.1, lwd=1) +
  geom_line(data = macro_output_mean, mapping = aes(x = time, y = mean_N), col="darkgreen", lwd = 1) +
	geom_line(data = macro_output_mean, mapping = aes(x = time, y = mean_P), col="darkred", lwd = 1) +
	geom_line(data = macro_out, mapping = aes(x = time, y = N), col="darkgreen", lwd=1, linetype = 2) +
	geom_line(data = macro_out, mapping = aes(x = time, y = P), col="darkred", lwd=1, linetype = 2) +
  geom_line(data = macro_proj_out, mapping = aes(x = time, y = N), col="darkgreen", lwd=1, linetype = 3) +
	geom_line(data = macro_proj_out, mapping = aes(x = time, y = P), col="darkred", lwd=1, linetype = 3)
#	xlim(0, 40) + ylim(50, 100)


ggplot() +
	geom_line(data = macro_output_stepproj, mapping = aes(x=time, y = N, col = as.factor(run), group=start), alpha = 0.5) +
#	geom_line(data = macro_output_stepfull, mapping = aes(x=time, y = P, group = paste(run, start)), col = "red", malpha = 1) +
	geom_line(data = macro_proj_out, mapping = aes(x = time, y = N), col="green", lwd=2) +
	geom_point(data = macro_proj_out, mapping = aes(x = time, y = N), size = 5) +
#  geom_line(data = macro_out, mapping = aes(x = time, y = P), col="red", lwd=1) +
	xlim(0, 10) + ylim(70, 100)

parms$micro_record = file("micro4.txt", open="w+")

macro_ab_out = macro_state_run_ab(macro_state, parms, 0)
close(parms$micro_record)

macro_ab_out = tbl_df(as.data.frame(macro_ab_out))
names(macro_ab_out) = c("time", "N", "P", "dN", "dP", "ddN", "ddP")

ggplot() +
  geom_line(data = macro_output_short, mapping = aes(x = time, y = N, group = run), col="green", alpha = 0.1, lwd=1) +
  geom_line(data = macro_output_short, mapping = aes(x = time, y = P, group = run), col="red", alpha = 0.1, lwd=1) +
  geom_line(data = macro_output_mean, mapping = aes(x = time, y = mean_N), col="darkgreen", lwd = 1) +
	geom_line(data = macro_output_mean, mapping = aes(x = time, y = mean_P), col="darkred", lwd = 1) +
	geom_line(data = macro_out, mapping = aes(x = time, y = N), col="darkgreen", lwd=1, linetype = 2) +
	geom_line(data = macro_out, mapping = aes(x = time, y = P), col="darkred", lwd=1, linetype = 2) +
  geom_line(data = macro_proj_out, mapping = aes(x = time, y = N), col="darkgreen", lwd=1, linetype = 3) +
	geom_line(data = macro_proj_out, mapping = aes(x = time, y = P), col="darkred", lwd=1, linetype = 3) +
  geom_line(data = macro_ab_out, mapping = aes(x = time, y = N), col="darkgreen", lwd=1, linetype = 4) +
  geom_line(data = macro_ab_out, mapping = aes(x = time, y = P), col="darkred", lwd=1, linetype = 4)


ggplot() +
	geom_line(data = macro_output_short, mapping = aes(x = time, y = N, group = run), col="green", alpha = 0.1, lwd=1) +
#	geom_line(data = macro_output_short, mapping = aes(x = time, y = P, group = run), col="red", alpha = 0.1, lwd=1) +
  geom_line(data = macro_output_mean, mapping = aes(x = time, y = mean_N), col="darkgreen", lwd = 1) +
#	geom_line(data = macro_output_mean, mapping = aes(x = time, y = mean_P), col="darkred", lwd = 1) +
	geom_line(data = macro_out, mapping = aes(x = time, y = N), col="darkgreen", lwd=1, linetype = 2) +
#	geom_line(data = macro_out, mapping = aes(x = time, y = P), col="darkred", lwd=1, linetype = 2) +
  geom_line(data = macro_proj_out, mapping = aes(x = time, y = N), col="darkgreen", lwd=1, linetype = 3) +
#	geom_line(data = macro_proj_out, mapping = aes(x = time, y = P), col="darkred", lwd=1, linetype = 3) +
  geom_line(data = macro_ab_out, mapping = aes(x = time, y = N), col="darkgreen", lwd=1, linetype = 4)
	#geom_line(data = macro_ab_out, mapping = aes(x = time, y = P), col="darkred", lwd=1, linetype = 4)



micro_output_ab = load_micro_output('micro4.txt')
macro_output_ab = micro_output_ab %>%
  arrange(start, time, run) %>%
  group_by(run, start, time) %>%
  summarise(N = sum(population), P = sum(infections*population))


ggplot() +
	geom_line(data = macro_output_ab, mapping = aes(x=time, y = N, col = as.factor(run), group=paste(run, start)), alpha = 0.1, lwd=1) +
#	geom_line(data = macro_output_stepfull, mapping = aes(x=time, y = P, group = paste(run, start)), col = "red", malpha = 1) +
	geom_line(data = macro_ab_out, mapping = aes(x = time, y = N), col="green", lwd=2) +
	geom_point(data = macro_ab_out, mapping = aes(x = time, y = N), size = 5) +
#  geom_line(data = macro_out, mapping = aes(x = time, y = P), col="red", lwd=1) +
	xlim(0, 20) + ylim(70, 100)

parms$micro_record = file("micro5.txt", open="w+")
macro_rk_out = macro_state_run_rk(macro_state, parms, 0)
close(parms$micro_record)

macro_rk_out = tbl_df(as.data.frame(macro_rk_out))
names(macro_rk_out) = c("time", "N", "P", "dN", "dP", "ddN", "ddP")

micro_output_rk = load_micro_output('micro5.txt')
macro_output_rk = micro_output_rk %>%
  arrange(start, time, run) %>%
  group_by(run, start, time) %>%
  summarise(N = sum(population), P = sum(infections*population))

ggplot() +
  geom_line(data = macro_output_short, mapping = aes(x = time, y = N, group = run), col="green", alpha = 0.1, lwd=1) +
  geom_line(data = macro_output_short, mapping = aes(x = time, y = P, group = run), col="red", alpha = 0.1, lwd=1) +
  geom_line(data = macro_output_mean, mapping = aes(x = time, y = mean_N), col="darkgreen", lwd = 1) +
  geom_line(data = macro_output_mean, mapping = aes(x = time, y = mean_P), col="darkred", lwd = 1) +
	geom_line(data = macro_out, mapping = aes(x = time, y = N), col="darkgreen", lwd=1, linetype = 2) +
	geom_line(data = macro_out, mapping = aes(x = time, y = P), col="darkred", lwd=1, linetype = 2) +
  geom_line(data = macro_proj_out, mapping = aes(x = time, y = N), col="darkgreen", lwd=1, linetype = 3) +
	geom_line(data = macro_proj_out, mapping = aes(x = time, y = P), col="darkred", lwd=1, linetype = 3) +
  geom_line(data = macro_ab_out, mapping = aes(x = time, y = N), col="darkgreen", lwd=1, linetype = 4) +
  geom_line(data = macro_ab_out, mapping = aes(x = time, y = P), col="darkred", lwd=1, linetype = 4) +
  geom_line(data = macro_rk_out, mapping = aes(x = time, y = N), col="darkgreen", lwd=2, linetype = 1) +
  geom_line(data = macro_ab_out, mapping = aes(x = time, y = P), col="darkred", lwd=2, linetype = 1)


