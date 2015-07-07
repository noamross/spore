micro_state.step = function(micro_state, parms, time) {
	is = seq_along(micro_state) - 1
	infections = sum(is*micro_state)
  N = sum(micro_state)
	if (N == 0) return(list(micro_state = micro_state, time_next = Inf))
	rates = c(
		P_inf = (parms$lambda * infections + parms$lambda_ex)*N,
		P_rec = parms$mu * infections,
		P_die = parms$alpha * infections + (parms$d * N),
		P_birth = max(parms$r * N * (1 - N/parms$K), 0)
	)

	time_next = time + rexp(1, sum(rates))
	event = sample(1:4, 1, prob=rates)

		switch(event, {
		i = sample(is, 1, prob=micro_state)
		micro_state[i + 1] = micro_state[i + 1] - 1
		if (i < length(micro_state) - 1) {
			micro_state[i + 2] = micro_state[i + 2] + 1
		}
	}, {
		i = sample(is, 1, prob = parms$mu * is * micro_state)
		micro_state[i+1] = micro_state[i+1] - 1
		micro_state[i] = micro_state[i] + 1
	}, {
		i = sample(is, 1, prob=parms$alpha * is * micro_state + parms$d*micro_state)
		micro_state[i+1] = micro_state[i+1] - 1
	}, {
		micro_state[1] = micro_state[1] + 1
	})
	return(list(micro_state = micro_state, time_next = time_next))
}

rates.micro_state = function(micro_state, parms) {
	is = seq_along(micro_state) - 1
	infections = sum(is*micro_state)
  N = sum(micro_state)

	rates = c(
		P_inf = (parms$lambda * infections + parms$lambda_ex)*N,
		P_rec = parms$mu * infections,
		P_die = parms$alpha * infections + (parms$d * N),
		P_birth = parms$r * N * (1 - N/parms$K)
	)
	return(rates)
	}

micro_state.stepto = function(micro_state, parms, time, timeto, record=NULL, run = 1) {

	time0 = time

	if(!is.null(record)) {
		micro_state.record(micro_state, time, connection = record, run = run, start = time0)
	}


	while(time < timeto) {
		out = micro_state.step(micro_state, parms, time)
		micro_state = out$micro_state
		time = out$time_next

		if(!is.null(record)) {
			micro_state.record(micro_state, time, connection = parms$micro_record, run = run, start = time0)
		}

	}
	return(micro_state)
}

micro_state.record = function(micro_state, time, connection, run, start) {
	writeChar(paste(c(run, start, time, micro_state, "\n"), collapse=" "), con=connection, eos=NULL)
}


lift.macro_state_FULL = function(macro_state, parms) {
  macro_state = pmax.int(macro_state, c(0,0))
	vals = rpois(floor(macro_state[1]) + (runif(1) < (macro_state[1] - floor(macro_state[1]))), macro_state[2]/macro_state[1])
	tabulate(bin = vals + 1, nbins = parms$max_i + 1)
}

lift.macro_state = cmpfun(lift.macro_state_FULL, options = list(optimize = 3))

restrict.micro_state = function(micro_state) {
 c(sum(micro_state), 	sum(micro_state * (seq_along(micro_state) -1)))
}
