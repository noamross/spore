Rprof("out_ad.out")

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
	t_max = 300
)

t = 0
pop = rep(0, parms$max_i + 1)
is = 0:(parms$max_i)
pop[1] = parms$init_pop

fname = "out_ad.txt"
if(file.exists(fname)) file.remove(fname)
outfile = file(fname, "w")
writeChar(paste(c(t, pop, "\n"), collapse=" "), con=outfile, eos=NULL)

max_rates_per_capita = c(
	P_inf_max = parms$lambda_ex + parms$lambda * 200, #parms$max_i * parms$K  + ,
	P_rec_max = parms$mu * parms$max_i,
	P_die_max = parms$alpha * parms$max_i + parms$d,
	P_birth_max = parms$r)
sum_max_rates_per_capita = sum(max_rates_per_capita)

infections = sum(is*pop)
N = sum(pop)
t_next = t
while (t < parms$t_max) {
	if (N == 0) break()

 	t_next = t_next + rexp(1, sum_max_rates_per_capita*N)
	i = sample(is, 1, prob=pop)
	event = sample(1:4, 1, prob=max_rates_per_capita)
	occurs = runif(1) < switch(event, 
		(parms$lambda * infections + parms$lambda_ex)/max_rates_per_capita[1],
		(parms$mu * (i - 1))/max_rates_per_capita[2],
		(parms$alpha * (i - 1) + parms$d)/max_rates_per_capita[3],
		(parms$r * (1 - N/parms$K))/max_rates_per_capita[4])
		
	if(occurs) {
	switch(event, {
		if (i < length(pop) - 1) {
			pop[i + 1] = pop[i + 1] - 1
			pop[i + 2] = pop[i + 2] + 1
		}
	}, {
		pop[i+1] = pop[i+1] - 1
		pop[i] = pop[i] + 1
	}, {
		pop[i+1] = pop[i+1] - 1
	}, {
		pop[1] = pop[1] + 1
	})
	
	t = t_next
	writeChar(paste(c(t, pop, "\n"), collapse=" "), con=outfile, eos=NULL)
	infections = sum(is*pop)
	N = sum(pop)

	}
	#browser()
}
close(outfile)
Rprof(NULL)

