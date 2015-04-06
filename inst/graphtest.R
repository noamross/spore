library(dplyr)
library(magrittr)

parms = list(
	max_i = 20,
	lambda = 0.001,
	lambda_ex = 0.1,
	alpha = 0.01,
	mu = 0.01,
	r = 0.5,
	d = 0.1,
	K = 100,
	init_pop = 100,
	time_max = 75,
	prevent_inf = 0,
	prevent_ex = 0,
	macro_timestep = 2,
	micro_timestep = 0.1,
	micro_relax_steps = 3,
	project = FALSE,
	n_sims = 50,
	micro_record = file("micro.txt", open="w+")
	#	macro_record = file("macro.txt", open="w")
)

micro_state = c(100, rep(0, parms$max_i))
macro_state = restrict.micro_state(micro_state)

for(run in 1:50) {
	micro_state.stepto(micro_state, parms = parms, time = 0, timeto= 75, record = parms$micro_record, run = run)
}
close(parms$micro_record)

micro_output = load_micro_output('micro.txt')

macro_output_short = restrict_micro_output(micro_output) %>%
	group_by(run) %>%
	mutate(time = plyr::round_any(time, 1, f=floor)) %>%
	filter(!duplicated(time)) %>%
	group_by()

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
close(parms$micro_record)
names(macro_out) = c("time", "N", "P")


micro_output_stepfull = load_micro_output('micro2.txt') 
macro_output_stepfull = restrict_micro_output(micro_output_stepfull)
macro_output_stepfull %<>% arrange(start, time, run)

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
macro_output_stepproj = restrict_micro_output(micro_output_stepproj)
macro_output_stepproj %<>% arrange(start, time, run)

ggplot() +
	geom_line(data = macro_output_short, mapping = aes(x = time, y = N, group = run), col="green", alpha = 0.1, lwd=1) +
#	geom_line(data = macro_output_short, mapping = aes(x = time, y = P, group = run), col="red", alpha = 0.1, lwd=1) +
	geom_line(data = macro_output_mean, mapping = aes(x = time, y = mean_N), col="green", lwd = 1) +
#	geom_line(data = macro_output_mean, mapping = aes(x = time, y = mean_P), col="red", lwd = 1) +
	geom_line(data = macro_out, mapping = aes(x = time, y = N), col="green", lwd=1, linetype = 2) +
#	geom_line(data = macro_out, mapping = aes(x = time, y = P), col="red", lwd=1, linetype = 2) +
  geom_line(data = macro_proj_out, mapping = aes(x = time, y = N), col="green", lwd=1, linetype = 3) +
#	geom_line(data = macro_proj_out, mapping = aes(x = time, y = P), col="red", lwd=1, linetype = 3) +
	xlim(0, 40) + ylim(50, 100)


ggplot() +
	geom_line(data = macro_output_stepproj, mapping = aes(x=time, y = N, col = as.factor(run), alpha = 0.5)) +
#	geom_line(data = macro_output_stepfull, mapping = aes(x=time, y = P, group = paste(run, start)), col = "red", alpha = 1) +
	geom_line(data = macro_proj_out, mapping = aes(x = time, y = N), col="green", lwd=1) +
	geom_point(data = macro_proj_out, mapping = aes(x = time, y = N)) +
#  geom_line(data = macro_out, mapping = aes(x = time, y = P), col="red", lwd=1) +
	xlim(0, 2) + ylim(70, 100)

