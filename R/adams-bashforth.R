macro_state_run_ab = function(macro_state, parms, time) {

  times = seq(time, parms$time_max, parms$macro_timestep)
  alpha = 1 + (parms$macro_timestep - parms$micro_timestep*parms$micro_relax_steps)/(2*parms$macro_timestep)
  project_timestep = parms$macro_timestep - parms$micro_timestep*(1 + parms$micro_relax_steps)

  output = cbind(times, matrix(NA, ncol = length(macro_state)*3, nrow=length(times)))
  output[1, -1] = c(macro_state, rep(NA, length(macro_state)*2))
  for(step in seq_along(times)) {

    time = times[step]

    if(step == 1) {
      lastvals = sapply(1:parms$n_sims, function(run) {
        micro_state = lift.macro_state(macro_state)
        relaxed_time = time + parms$micro_timestep*parms$micro_relax_steps
        micro_state_relaxed = micro_state.stepto(micro_state, parms, time = time, timeto = relaxed_time, run = run, record=parms$micro_record)
        next_time = relaxed_time + parms$micro_timestep
        micro_state_next = micro_state.stepto(micro_state_relaxed, parms, time = relaxed_time, timeto = next_time, run = run, record=parms$micro_record)
        micro_state_full = micro_state.stepto(micro_state_next, parms, time = next_time, timeto = times[2], run = run, record=parms$micro_record)
        macro_state_relaxed = restrict.micro_state(micro_state_relaxed)
        macro_state_next = restrict.micro_state(micro_state_next)
        macro_state_full = restrict.micro_state(micro_state_full)
        return(c(macro_state_relaxed, macro_state_next, macro_state_full))
      })
      last_deriv = (rowMeans(lastvals[3:4,]) - rowMeans(lastvals[1:2,]))/parms$micro_timestep
      macro_state = rowMeans(lastvals[5:6,])
      macro_deriv = last_deriv
      macro_second_deriv = rep(NA, length(macro_state))
    } else {
      vals = sapply(1:parms$n_sims, function(run) {
        micro_state = lift.macro_state(macro_state)
        relaxed_time = time + parms$micro_timestep*parms$micro_relax_steps
        micro_state_relaxed = micro_state.stepto(micro_state, parms, time = time, timeto = relaxed_time, run = run, record=parms$micro_record)
        next_time = relaxed_time + parms$micro_timestep
        micro_state_next = micro_state.stepto(micro_state_relaxed, parms, time = relaxed_time, timeto = next_time, run = run, record=parms$micro_record)
        macro_state_relaxed = restrict.micro_state(micro_state_relaxed)
        macro_state_next = restrict.micro_state(micro_state_next)
        return(c(macro_state_relaxed, macro_state_next))
      })
      macro_state_next = rowMeans(vals[3:4,])
      next_deriv = (macro_state_next - rowMeans(vals[1:2,]))/parms$micro_timestep
      macro_state = macro_state_next + (alpha*next_deriv + (1 - alpha)*last_deriv)*project_timestep
      macro_deriv = (project_timestep*next_deriv + (parms$micro_relax_steps*parms$micro_timstep + 1)*last_deriv)/parms$macro_timestep
      macro_second_deriv = (next_deriv - last_deriv)/parms$macro_timestep
      last_deriv = next_deriv
    }
    output[step, -seq(1:(1+length(macro_state)))] = c(macro_deriv, macro_second_deriv)
    if(step != length(times)) {
      output[step + 1, 2:(1+length(macro_state))] = macro_state
    }
  }
  return(output)
}




