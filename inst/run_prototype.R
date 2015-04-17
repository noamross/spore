 Rprof("out_gil.out")

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
pop = rep(0, parms$max_i + 1)
is = 0:(parms$max_i)
pop[1] = parms$init_pop

fname = "out_gil.txt"
if(file.exists(fname)) file.remove(fname)
outfile = file(fname, "w")
writeChar(paste(c(t, pop, "\n"), collapse=" "), con=outfile, eos=NULL)

while (t < parms$t_max) {
	infections = sum(is*pop)
  N = sum(pop)
	if (N == 0) break()
	rates = c(
		P_inf = (parms$lambda * infections + parms$lambda_ex)*N,
		P_rec = parms$mu * infections,
		P_die = parms$alpha * infections + (parms$d * N),
		P_birth = parms$r * N * (1 - N/parms$K)
	)
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

	#browser()
	t = t_next
	writeChar(paste(c(t, pop, "\n"), collapse=" "), con=outfile, eos=NULL)
}
close(outfile)
Rprof(NULL)
