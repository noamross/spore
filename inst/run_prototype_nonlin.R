parms = list(
	max_i = 200,
	lambda = 0,
	lambda_ex = 0.1,
	alpha = 0.01,
	mu = 0,
	r = 0,
	d = 0,
	K = 1000,
	init_pop = 1000,
	t_max = 1000
)

t = 0
pop = rep(0, parms$max_i + 1)
is = 0:(parms$max_i)
pop[1] = parms$init_pop

fname = "out.txt"
cat(c(t, pop, "\n"), file=fname, append=FALSE)

while (t < parms$t_max) {
	infections = sum(is*pop)
  N = sum(pop)
	if (N == 0) break()
	rates = c(
		P_inf = (parms$lambda * infections + parms$lambda_ex)*N,
		P_rec = parms$mu * infections,
		P_die = parms$alpha * sum(is * pop) + (parms$d * N),
		P_birth = parms$r * N * (1 - N/parms$K)
	)
	t_next = t + rexp(1, sum(rates))
  event = sample(1:4, 1, prob=rates)	
  if (event == 4) {
  	pop[1] = pop[1] + 1
  } else if (event == 3) {
		i = sample(is, 1, prob=parms$alpha * (is * pop) + parms$d*pop)
		pop[i+1] = pop[i+1] - 1
  } else if (event == 2) {
  	i = sample(is, 1, prob = parms$mu * is * pop)
  	pop[i+1] = pop[i+1] - 1
  	pop[i] = pop[i] + 1
  } else if (event == 1) {
  	i = sample(is, 1, prob=pop)
  	if (i == length(pop) - 1) {
  		pop[i + 2] = 0
  		is[i + 2] = i + 1
  	}
  	pop[i + 1] = pop[i + 1] - 1
  	pop[i + 2] = pop[i + 2] + 1
  }
	#browser()
	cat(c(t_next, pop, "\n"), file=fname, append=TRUE)
	t = t_next
}
