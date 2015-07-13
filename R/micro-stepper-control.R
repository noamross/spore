#' @export
micro_state_c.step_FULL = function(micro_state, parms, control, time) {
	is = seq_along(micro_state) - 1
	infections = sum(is*micro_state)
  N = sum(micro_state)
	if (N == 0) return(list(micro_state = micro_state, time_next = Inf))
	rates = c(
		P_inf = (parms$lambda * infections +
		           parms$lambda_ex * exp(-control))*N,
		P_rec = parms$mu * infections,
		P_die = parms$alpha * (infections^parms$alpha_power) + (parms$d * N),
		P_birth = max(parms$r * N * (1 - N/parms$K), 0)
	)

	time_next = time + rexp(1, sum(rates))
	event = sample.int(4, 1, prob=rates)

		switch(event, {
		i = sample.int(parms$max_i + 1, 1, prob=micro_state) - 1
		micro_state[i + 1] = micro_state[i + 1] - 1
		if (i < length(micro_state) - 1) {
			micro_state[i + 2] = micro_state[i + 2] + 1
		}
	}, {
		i = sample.int(parms$max_i + 1, 1, prob = parms$mu * is * micro_state) - 1
		micro_state[i+1] = micro_state[i+1] - 1
		micro_state[i] = micro_state[i] + 1
	}, {
		i = sample.int(parms$max_i + 1, 1, prob=parms$alpha * (is^parms$alpha_power) * micro_state + parms$d*micro_state) - 1
		micro_state[i+1] = micro_state[i+1] - 1
	}, {
		micro_state[1] = micro_state[1] + 1
	})
	return(list(micro_state = micro_state, time_next = time_next))
}



#' @importFrom compiler cmpfun
#' @export
micro_state_c.step = cmpfun(micro_state_c.step_FULL, options = list(optimize = 3))

rates.micro_state = function(micro_state, parms, control) {
	is = seq_along(micro_state) - 1
	infections = sum(is*micro_state)
  N = sum(micro_state)

	rates = c(
		P_inf = (parms$lambda * infections +
		           parms$lambda_ex * exp(-control))*N,
		P_rec = parms$mu * infections,
		P_die = parms$alpha * infections^parms$alpha_power + (parms$d * N),
		P_birth = parms$r * N * (1 - N/parms$K)
	)
	return(rates)
	}

#' @export
micro_state_c.stepto = function(micro_state, parms, control, time, timeto, record=NULL, run = 1) {

	time0 = time
  RECORD = !is.null(record) && ("connection" %in% class(parms$micro_record))

	if(RECORD) {
		micro_state_c.record(micro_state, control, time, connection = record, run = run, start = time0)
	}


	while(time < timeto) {
		out = micro_state_c_step(micro_state, parms, control, time)
		if(out$time_next < timeto) micro_state = out$micro_state
		time = out$time_next

		if(RECORD) {
			micro_state_c.record(micro_state, control, time, connection = record, run = run, start = time0)
		}

	}
	return(micro_state)
}

micro_state_cpath.stepto = function(micro_state, parms, controlfn, time, timeto, record=NULL, run = 1) {

  time0 = time
  control = controlfn(time)
  if(!is.null(record)) {
    micro_state_c.record(micro_state, control, time, connection = record, run = run, start = time0)
  }


  while(time < timeto) {
    control = controlfn(time)
    out = micro_state_c_step(micro_state, parms, control, time)
    micro_state = out$micro_state
    time = out$time_next

    if("connection" %in% class(parms$micro_record)) {
      micro_state_c.record(micro_state, control, time, connection = parms$micro_record, run = run, start = time0)
    }

  }
  return(micro_state)
}

micro_state_c.record = function(micro_state, control, time, connection, run, start) {
	writeChar(paste(c(run, start, time, control, micro_state, "\n"), collapse=" "), con=connection, eos=NULL)
}

