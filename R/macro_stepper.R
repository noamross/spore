macro_state.stepfull = function(macro_state, parms, time) {
		macro_states = sapply(1:parms$n_sims, function(run) {
			micro_state = lift.macro_state(macro_state)
      micro_state = micro_state.stepto(micro_state, parms, time = time, timeto = time + parms$macro_timestep, record = parms$micro_record, run = run)
			return(restrict.micro_state(micro_state))
		})
		macro_state = rowMeans(macro_states)
		return(macro_state)
}

macro_state.stepproject = function(macro_state, parms, time) {
		macro_states = sapply(1:parms$n_sims, function(run) {
			micro_state = lift.macro_state(macro_state)
		  relaxed_time = time + parms$micro_timestep*parms$micro_relax_steps
      micro_state_relaxed = micro_state.stepto(micro_state, parms, time = time, timeto = relaxed_time, run = run)
      next_time = relaxed_time + parms$micro_timestep
      micro_state_next = micro_state.stepto(micro_state_relaxed, parms, time = relaxed_time, timeto = next_time, run = run)
      
			macro_state_relaxed = restrict.micro_state(micro_state_relaxed)
			macro_state_next = restrict.micro_state(micro_state_next)
			
			return(macro_state_relaxed + ((macro_state_next - macro_state_relaxed) /
																		(next_time - relaxed_time)) *
						 	                     (time + parms$macro_timestep - next_time))
		})
		macro_state = rowMeans(macro_states)
		return(macro_state)
}

macro_state.run = function(macro_state, parms, time) {
	
	times = seq(time, parms$time_max, parms$macro_timestep)
	
	output = cbind(times, matrix(NA, ncol = length(macro_state), nrow=length(times)))
	output[1, -1] = macro_state
	for(step in seq_along(times)[-1]) {
	
	time = times[step - 1]
	
	if(parms$project) {
		macro_state = macro_state.stepproject(macro_state, parms, time)
	} else {
		macro_state = macro_state.stepfull(macro_state, parms, time)
	}
	
	output[step, -1] = macro_state

	}
	return(output)
}
