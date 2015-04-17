library(polynom)

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
	t_max = 300,
	prevent_inf = 0,
	prevent_ex = 0,
	micro_timestep = 0.1,
	micro_relax_steps = 3,
	project = TRUE,
	n_sims = 10,
	out_file = "outfile.txt"
)

t = 0
t_macro = 0
macro_pop = c(parms$init_pop, 0)
pop = rep(0, parms$max_i + 1)
is = 0:(parms$max_i)
pop[1] = parms$init_pop
macro_times = seq(0, parms$t_max, parms$macro_timestep)
macro_step = 1
# fname = "out_gil.txt"
# if(file.exists(fname)) file.remove(fname)
# outfile = file(fname, "w")
# writeChar(paste(c(t, pop, "\n"), collapse=" "), con=outfile, eos=NULL)
macro_vals = rbind(macro_times, macro_times)
macro_vals[] = NA

round_rand = function(x) {
  prob = x - floor(x)
  return(ifelse(runif(length(x)) < prob, ceiling(x), floor(x)))
	}

lift = function(macro_pop) {
	vals = rpois(round_rand(macro_pop[1]), macro_pop[2]/macro_pop[1])
	as.integer(table(factor(vals, levels = 0:max(c(parms$max_i, max(vals))))))
}

restrict = function(pop) {
 c(sum(pop), 	sum(pop * (seq_along(pop) -1)))
}

change_indices = function(x) {
	ind = which(x[2:length(x)] - x[1:(length(x) - 1)] != 0) + 1
	return(c(ind[1] - 1, ind))
	}

for (macro_step in seq_along(macro_times)) {
	
	t_macro = macro_times[macro_step]
	sims_N = rep(NA, parms$n_sims)
	sims_P = rep(NA, parms$n_sims)
	for(sim in 1:parms$n_sims) {
		
		k_N = rep(NA, 2*parms$micro_steps_meas)
		k_P = rep(NA, 2*parms$micro_steps_meas)
		k_times = rep(NA, 2*parms$micro_steps_meas)
		
		
		pop = lift(macro_pop)
		t = t_macro
		is = seq_along(pop) - 1
		infections = sum(is*pop)
		N = sum(pop)
		
		relaxed = FALSE
		micro_step = 0
		while(!relaxed) {
			micro_step = micro_step + 1
			rates = c(
				P_inf = (parms$lambda * infections * exp(-parms$prevent_inf) + parms$lambda_ex * exp(-parms$prevent_ex))*N,
				P_rec = parms$mu * infections,
				P_die = parms$alpha * infections + (parms$d * N),
				P_birth = max(0, parms$r * N * (1 - N/parms$K)))
			
			t_next = t + rexp(1, sum(rates))
			event = sample(1:4, 1, prob=rates)	
			
			switch(event, {
				i = sample(is, 1, prob=pop)
				if (i < length(pop) - 1) {
					pop[i + 1] = pop[i + 1] - 1
					pop[i + 2] = pop[i + 2] + 1
				}
			}, {
				i = sample(is, 1, prob = parms$mu * is * pop)
				pop[i+1] = pop[i+1] - 1
				pop[i] = pop[i] + 1
			}, {
				i = sample(is, 1, prob=parms$alpha * is * pop + parms$d*pop)
				pop[i+1] = pop[i+1] - 1
			}, {
				pop[1] = pop[1] + 1
			})
			
			macro_state_tmp = restrict(pop)
			k_N[micro_step] = macro_state_tmp[1]
			k_P[micro_step] = macro_state_tmp[2]
			k_times[micro_step] = t_next
			t = t_next
			if(length(change_indices(k_N)) > parms$micro_steps_init + parms$micro_steps_meas &
				 length(change_indices(k_P)) > parms$micro_steps_init + parms$micro_steps_meas &
				 micro_step > parms$micro_steps_meas) {
				relaxed = TRUE
			}
			#	writeChar(paste(c(sim, t, pop, "\n"), collapse=" "), con=outfile, eos=NULL)
		}
		#print(rbind(k_N, k_P, k_times))
		k_N = na.omit(k_N)
		k_P = na.omit(k_P)
		k_times = na.omit(k_times)
		
		
	  sims_N[sim] = k_N[micro_step] + 
			((k_N[micro_step] - k_N[change_indices(k_N)[parms$micro_steps_init + 1]]) / 
	  	 (k_times[micro_step] - k_times[change_indices(k_N)[parms$micro_steps_init + 1]])) *
			   (macro_times[macro_step + 1] - k_times[micro_step])
	  sims_P[sim] = k_P[micro_step] + 
			((k_P[micro_step] - k_P[change_indices(k_P)[parms$micro_steps_init + 1]]) / 
	  	 (k_times[micro_step] - k_times[change_indices(k_P)[parms$micro_steps_init + 1]])) *
			   (macro_times[macro_step + 1] - k_times[micro_step])
	}
	macro_pop = pmax(c(mean(sims_N), mean(sims_P)), c(0,0))
	macro_vals[,macro_step] = macro_pop
}

close(outfile)
Rprof(NULL)
	
library(polynom)
XX = 1:3
YY = c(1,2,1)
plot(XX,YY)
p = poly.calc(x= XX, y=YY)
pf = as.function(p)
plot(p)
predict(p, 0.5)
points(c(XX, 0.5), c(YY, p(0.5)))
