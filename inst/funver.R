micro_state.step = function(micro_state, parms, time) {
	is = seq_along(micro_state) - 1
	infections = sum(is*micro_state)
  N = sum(micro_state)
	if (N == 0) return(list(micro_state = micro_state, time_next = Inf))
	rates = c(
		P_inf = (parms$lambda * infections + parms$lambda_ex)*N,
		P_rec = parms$mu * infections,
		P_die = parms$alpha * infections + (parms$d * N),
		P_birth = parms$r * N * (1 - N/parms$K)
	)
#browser()
	time_next = time + rexp(1, sum(rates))
	event = sample(1:4, 1, prob=rates)	
#browser()	
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

micro_state.stepto = function(micro_state, parms, time, timeto, record=NULL) {
	if(!is.null(record)) {
		micro_state.record(micro_state, time, connection = record)
	}
	
	while(time < timeto) {
		out = micro_state.step(micro_state, parms, time)
		micro_state = out$micro_state
		time = out$time_next
		
		if(!is.null(record)) {
			micro_state.record(micro_state, time, connection = record)
		}
		
	}
}

micro_state.record = function(micro_state, time, connection) {
	writeChar(paste(c(time, micro_state, "\n"), collapse=" "), con=connection, eos=NULL)
}
