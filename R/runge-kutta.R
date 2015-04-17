macro_state_run_rk = function(macro_state, parms, time) {

  times = seq(time, parms$time_max, parms$macro_timestep)
  alpha = 1/2 + (parms$micro_timestep*parms$micro_relax_steps)/(2*parms$macro_timestep)
  project_timestep = parms$macro_timestep - parms$micro_timestep*(1 + parms$micro_relax_steps)

  output = cbind(times, matrix(NA, ncol = length(macro_state)*3, nrow=length(times)))
  output[1, -1] = c(macro_state, rep(NA, length(macro_state)*2))
  for(step in seq_along(times)) {

    time = times[step]

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
      deriv = (macro_state_next - rowMeans(vals[1:2,]))/parms$micro_timestep
      est_forward_macro_state = macro_state_next + deriv*project_timestep
      forward_time = time + parms$macro_timestep
    nextvals = sapply(1:parms$n_sims, function(run) {
      micro_state = lift.macro_state(est_forward_macro_state)
      relaxed_time = forward_time + parms$micro_timestep*parms$micro_relax_steps
      micro_state_relaxed = micro_state.stepto(micro_state, parms, time = forward_time, timeto = relaxed_time, run = run, record=NULL)
      next_time = relaxed_time + parms$micro_timestep
      micro_state_next = micro_state.stepto(micro_state_relaxed, parms, time = relaxed_time, timeto = next_time, run = run, record=NULL)
      macro_state_relaxed = restrict.micro_state(micro_state_relaxed)
      macro_state_next = restrict.micro_state(micro_state_next)
      return(c(macro_state_relaxed, macro_state_next))
    })
      macro_state_forward_next = rowMeans(nextvals[3:4,])
      forward_deriv = (macro_state_forward_next - rowMeans(nextvals[1:2,]))/parms$micro_timestep
      macro_state_projected = macro_state_next + (alpha*deriv + (1 - alpha)*forward_deriv)*project_timestep

    if(step != 1) {
      macro_deriv = (project_timestep*deriv + (parms$micro_relax_steps*parms$micro_timestep)*last_deriv)/parms$macro_timestep
      macro_second_deriv = (deriv - last_deriv)/parms$macro_timestep
      output[step, -seq(1:(1+length(macro_state)))] = c(macro_deriv, macro_second_deriv)
    }

    macro_state = macro_state_projected
    last_deriv = deriv

    if(step != length(times)) {
      output[step + 1, 2:(1+length(macro_state))] = macro_state_projected
    }
  }
  return(output)
}




