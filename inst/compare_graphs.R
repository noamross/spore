## ----ODEsim------------------------------------------------------------

library(deSolve)
syseq = function(t, state, parms) {
  with(as.list(c(state, parms)), {
    h_calc = suppressWarnings(log(- S2 * N * lambda_ex/c))
    if(is.nan(h_calc)) h = 0 else h = h_calc
    h = min(max(h, control_min), control_max)
    H_calc = H(h, t, state, parms)
    dN = max(r*N*(1 - N/K), 0) - alpha*P - d*N
    if(isTRUE(all.equal(dN, 0))) dN = 0
    dP = P*(lambda*N - mu - d - alpha - alpha*(P)/N) + exp(-h)*N*lambda_ex
    dS1 = -v - S1*(r - d - 2*(r/K)*N) - S2*(lambda*P + alpha*(P^2/N^2) + exp(-h) * lambda_ex) + delta*S1
    dS2 = S1*alpha - S2*(lambda*N - mu - d - alpha - 2*alpha*P/N) + delta*S2
    derivs = c(dN=dN, dP=dP, dS1=dS1, dS2=dS2)
    return(list(derivs, c(derivs, h_calc=h_calc, H_calc = H_calc, h=h)))
  })
}


H = function(h, t, state, parms) {
  with(as.list(c(state, parms)), {
    val = v*N - c*h + S1*(r*N*(1 - N/K) - alpha*P - d*N) + S2*(lambda*P*N - mu*P - d*P - alpha*P - alpha*(P^2)/N + exp(-h)*N*lambda_ex)
    return(val)
  })
}

N0 = macro_state[1]
P0 = macro_state[2]
init = c(N = N0, P = P0, S1 = 0, S2 = 0)
times = seq(0,40, by=1)

parms0 = parms
parms0$control_max = 0

out0 = lsoda(init, times, syseq, parms0)

out0 = as.data.frame(out0)

## ----IBMsim-------------------------------------------------------------

parms$micro_record = file("micro.txt", open="w+")
for(run in 1:parms$n_comp_sims) {
	micro_state_c.stepto(micro_state, parms = parms, time = 0, timeto=parms$time_max, record = parms$micro_record, run = run, control = 0)
  if(interactive()) cat(run, "\r")
}
close(parms$micro_record)

output = read_delim("micro.txt", " ", col_names = FALSE)
colnames(output) = c("run", "start", "time", "control", as.character(0:(ncol(output)-5)))

micro_output_short = output %>%
  group_by(run) %>%
  filter((!duplicated(floor(time), fromLast=TRUE)) | time == 0) %>%
  group_by() %>%
  mutate(time = ceiling(time)) %>%
  mutate_each(funs(na_to_0), -run, -start, -time, -control) %>%
  filter(time <= 41)

macro_output_short = micro_output_short %>%
  arrange(run, time) %>%
  gather("infections", "population", -run, -start, -time, -control) %>%
  arrange(run, time, infections) %>%
  mutate(infections = as.integer(as.character(infections))) %>%
  group_by(run, time) %>%
  summarize(N = sum(population), P = sum(population*infections)) %>%
  group_by()

macro_output_mean = macro_output_short %>%
  group_by(time) %>%
	summarize(mean_N = mean(N), mean_P = mean(P), lower_N = mean_N - 2*sd(N),
            upper_N = mean_N + 2*sd(N), lower_P = mean_P - 2*sd(P),
            upper_P = mean_P + 2*sd(P))

macro_output_mean_no_extinctions = macro_output_short %>%
  group_by(run) %>%
  filter(!(N == 0)) %>%
  group_by(time) %>%
	summarize(mean_N = mean(N), mean_P = mean(P))


## ----EFsim-------------------------------------------------------------

shadow_state = c(0, 0)
time = 0
parms$control_max = 0
zz <- file.remove('out.txt')
zz <- file.create('out.txt')
parmvec = unlist(as.relistable(within(parms, {n_sims = 100000; n_sims_jacob = 100000; parallel_cores = 3; progress=interactive()})))
ef_sims  <- ode(y = init, times = 0:40, parms=parmvec, func = opt_derivs, method = "euler")
ef_sims2 = read.csv("out.txt", header=FALSE)
names(ef_sims2) = c("time", "N", "P", "S1", "S2", "dN", "dP", "dS1", "dS2", "ddNdN", "ddPdN", "ddNdP", "ddPdP", "h", "H")

## ----fig1------------------------------------------------------------

ODEplot = ggplot() +
  geom_line(data = out0, mapping = aes(x = time, y = N), col="blue", lwd=1, lty=1) +
  geom_line(data = out0, mapping = aes(x = time, y = P), col="red", lwd=1, lty=1) +
  annotate("text", x=c(5,23), y=c(600, 900), label = c("pathogen", "host"),
           color = c("red", "blue"), size=5) +
  annotate("text", x = 35, y = 10, label = "ODE model", size = 5) +
  labs(list(Title = "", x = "", y = "")) + theme_nr

IBMplot = ggplot() +
	geom_line(data = macro_output_short, mapping = aes(x = time, y = N, group = run), col="steelblue", alpha = 0.05, lwd=1) +
	geom_line(data = macro_output_short, mapping = aes(x = time, y = P, group = run), col="tomato", alpha = 0.05, lwd=1) +
	geom_line(data = macro_output_mean, mapping = aes(x = time, y = mean_N), col="steelblue", lwd = 1, lty=2) +
	geom_line(data = macro_output_mean, mapping = aes(x = time, y = mean_P), col="tomato", lwd = 1, lty=2) +
  geom_line(data = macro_output_mean, mapping = aes(x = time, y = lower_N), col="steelblue", lwd = 1, lty=3) +
	geom_line(data = macro_output_mean, mapping = aes(x = time, y = lower_P), col="tomato", lwd = 1, lty=3) +
  geom_line(data = macro_output_mean, mapping = aes(x = time, y = upper_N), col="steelblue", lwd = 1, lty=3) +
	geom_line(data = macro_output_mean, mapping = aes(x = time, y = upper_P), col="tomato", lwd = 1, lty=3) +
  annotate("text", x = 35, y = 10, label = "IBM model", size = 5) +
  labs(list(Title = "", x = "", y = "")) + theme_nr

EFplot = ggplot() +
	geom_line(data = as.data.frame(ef_sims), mapping = aes(x = times, y = N), col="darkblue", lwd=1, lty=6) +
  geom_line(data = as.data.frame(ef_sims), mapping = aes(x = times, y = P), col="darkred", lwd=1, lty=6) +
  annotate("text", x = 35, y = 10, label = "EF model", size = 5) +
  labs(list(Title = "", x = "", y = "")) + theme_nr


ALLplot = ggplot() +
  geom_line(data = out0, mapping = aes(x = time, y = N), col="blue", lwd=1, lty=1) +
  geom_line(data = out0, mapping = aes(x = time, y = P), col="red", lwd=1, lty=1) +
  geom_line(data = macro_output_mean, mapping = aes(x = time, y = mean_N), col="steelblue", lwd = 1, lty=2) +
	geom_line(data = macro_output_mean, mapping = aes(x = time, y = mean_P), col="tomato", lwd = 1, lty=2) +
  geom_line(data = as.data.frame(ef_sims), mapping = aes(x = times, y = N), col="darkblue", lwd=1, lty=6) +
  geom_line(data = as.data.frame(ef_sims), mapping = aes(x = times, y = P), col="darkred", lwd=1, lty=6) +
  annotate("text", x = 35, y = 10, label = "All", size = 5) +
  labs(list(Title = "", x = "", y = "")) + theme_nr


cowplot::plot_grid(ODEplot, IBMplot, EFplot, ALLplot, ncol=2,
                   labels = c('A', 'B', 'C', 'D'), hjust = -4, vjust = 1) +
  draw_label('Host / Population Size', angle = 90, fontfamily = 'Lato',
             size = 22, vjust = -23.6) +
  draw_label('Time', fontfamily = 'Lato', size = 22, vjust = 13.3)
