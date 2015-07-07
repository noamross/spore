macro_state_c.stepfull = function(macro_state, parms, control, time) {
  macro_states = sapply(1:parms$n_sims, function(run) {
    micro_state = lift.macro_state(macro_state, parms)
    micro_state = micro_state_c.stepto(micro_state, parms, control, time = time, timeto = time + parms$macro_timestep, record = parms$micro_record, run = run)
    return(restrict.micro_state(micro_state))
  })
  macro_state = rowMeans(macro_states)
  return(macro_state)
}

macro_state_c.stepproject = function(macro_state, parms, control, time) {
  macro_states = sapply(1:parms$n_sims, function(run) {
    micro_state = lift.macro_state(macro_state, parms)
    relaxed_time = time + parms$micro_timestep*parms$micro_relax_steps
    micro_state_relaxed = micro_state.stepto(micro_state, parms, control, time = time, timeto = relaxed_time, run = run)
    next_time = relaxed_time + parms$micro_timestep
    micro_state_next = micro_state.stepto(micro_state_relaxed, parms, control, time = relaxed_time, timeto = next_time, run = run)

    macro_state_relaxed = restrict.micro_state(micro_state_relaxed)
    macro_state_next = restrict.micro_state(micro_state_next)

    return(macro_state_relaxed + ((macro_state_next - macro_state_relaxed) /
                                    (next_time - relaxed_time)) *
             (time + parms$macro_timestep - next_time))
  })
  macro_state = rowMeans(macro_states)
  return(macro_state)
}



macro_state_c_runopt = function(macro_state_init, parms, shadow_state_init, time, control_guess_init) {


  alpha = 1 + (parms$macro_timestep - parms$micro_timestep*(parms$micro_relax_steps - 1))/(2*parms$macro_timestep)
  project_timestep = parms$macro_timestep - parms$micro_timestep*(parms$micro_relax_steps)

  times = seq(time, parms$time_max, parms$macro_timestep)

  macro_states = matrix(NA, length(times), length(macro_state_init))
  macro_derivs = macro_states
  macro_second_derivs = macro_states
  shadow_states = matrix(NA, length(times), length(shadow_state_init))
  shadow_derivs = shadow_states
  controls = matrix(NA, length(times), length(control_guess_init))
  hamiltonian = rep(NA, length(times))
  alt = rep(FALSE, length(times))
  macro_states[1,] = macro_state_init
  shadow_states[1,] = shadow_state_init

  step = 1
  if(parms$progress) p <- progress_estimated(length(times))
  time = time[step]

  opt = determine_control(macro_state = macro_state_init, parms = parms,
                          shadow_state = shadow_state_init, time = time,
                          control_guess = control_guess_init)
  controls[1,] = opt$control
  macro_derivs[1,] = opt$macro_state_deriv
  hamiltonian[1] = opt$hamiltonian
  macro_second_derivs[1,] = second_deriv_from_3pts(
    c(time, time + parms$micro_relax_steps*parms$micro_timestep, time + (parms$micro_relax_steps + 1) * parms$micro_timestep),
    rbind(macro_state_init, opt$macro_state_relaxed, opt$macro_state_next)
  )
  rbind(macro_state_init, opt$macro_state_relaxed, opt$macro_state_next)
  macro_states[2,] = opt$macro_state_relaxed + macro_derivs[1,] * project_timestep

  shadow_derivs[step, 1] = -(parms$v +
                               shadow_states[step, 1] * macro_second_derivs[step,1] / macro_derivs[step,1] +
                               shadow_states[step, 2] * macro_second_derivs[step,2] / macro_derivs[step,1])
  shadow_derivs[step,2] =
    -(shadow_states[step, 1] * macro_second_derivs[step ,1] / macro_derivs[step ,2] +
        shadow_states[step, 2] * macro_second_derivs[step ,2] / macro_derivs[step ,2])

  if(any(macro_derivs[step,] == 0)) {
    alt_shadow_derivs = alt_shadow_derivs_calc(parms=parms, macro_state = macro_states[step,],
                                               macro_deriv = macro_derivs[step,],
                                               last_deriv_est = NULL,
                                               macro_second_deriv = macro_second_derivs[step,],
                                               shadow_state = shadow_states[step,],
                                               control = controls[step,],
                                               diff_step = 1)

    shadow_derivs[step, macro_derivs[step,] == 0] = alt_shadow_derivs[macro_derivs[step,] == 0]
    alt[step] = TRUE
  }

  shadow_states[2,] = shadow_states[1,] + shadow_derivs[1,]*parms$macro_timestep

  last_deriv_est = opt$macro_state_deriv
  last_second_deriv_est = macro_second_derivs[1,]

  if(parms$progress) p$tick()$print()
  for(step in seq_along(times)[-1]) {

    time = times[step]

    opt = determine_control(macro_states[step,], parms, shadow_states[step,], time, controls[step - 1])

    controls[step,] = opt$control
    hamiltonian[step] = opt$hamiltonian
    macro_derivs[step,] = (project_timestep*opt$macro_state_deriv + (parms$micro_relax_steps*parms$micro_timestep)*last_deriv_est)/parms$macro_timestep
    local_second_deriv = second_deriv_from_3pts(
      c(time, time + parms$micro_relax_steps*parms$micro_timestep, time + (parms$micro_relax_steps + 1) * parms$micro_timestep),
      rbind(macro_states[step,], opt$macro_state_relaxed, opt$macro_state_next)
    )
    macro_second_derivs[step,] = (project_timestep*local_second_deriv + (parms$micro_relax_steps*parms$micro_timestep)*last_second_deriv_est)/parms$macro_timestep

    #macro_second_derivs[step,] = (opt$macro_state_deriv - last_deriv_est)/parms$macro_timestep

    shadow_derivs[step, 1] = -(parms$v +
                                 shadow_states[step, 1] * macro_second_derivs[step,1] / macro_derivs[step,1] +
                                 shadow_states[step, 2] * macro_second_derivs[step,2] / macro_derivs[step,1])
    shadow_derivs[step,2] =
      -(shadow_states[step, 1] * macro_second_derivs[step ,1] / macro_derivs[step ,2] +
          shadow_states[step, 2] * macro_second_derivs[step ,2] / macro_derivs[step ,2])

    if(any(macro_derivs[step,] == 0)) {
      alt_shadow_derivs = alt_shadow_derivs_calc(parms=parms, macro_state = macro_states[step,],
                                                 macro_deriv = macro_derivs[step,],
                                                 last_deriv_est = last_deriv_est,
                                                 macro_second_deriv = macro_second_derivs[step,],
                                                 shadow_state = shadow_states[step,],
                                                 control = controls[step,],
                                                 diff_step = 1)

      shadow_derivs[step, macro_derivs[step,] == 0] = alt_shadow_derivs[macro_derivs[step,] == 0]
      alt[step] = TRUE
    }

    if(step != length(times)) {
      macro_states[step + 1, ] = opt$macro_state_next + (alpha*opt$macro_state_deriv + (1 - alpha)*last_deriv_est)*project_timestep
      shadow_states[step + 1, ] = shadow_states[step,] + (0.5 * shadow_derivs[step - 1,] + 1.5 * shadow_derivs[step,]) * parms$macro_timestep
    }
    last_deriv_est2 = last_deriv_est
    last_deriv_est = opt$macro_state_deriv
    last_second_deriv_est = macro_second_derivs[step,]
    if(parms$progress) p$tick()$print()
  }
  return(cbind(times, macro_states, macro_derivs, macro_second_derivs, shadow_states, shadow_derivs, controls, hamiltonian, alt))
}

#' @import nloptr
#' @export
determine_control = function(macro_state, parms, shadow_state, time, control_guess) {

  Hamiltonian = function(control, macro_state, parms, shadow_state, time) {
    if (control < parms$control_min | control > parms$control_max) return(Inf)
    vals = sapply(X = 1:parms$n_sims, FUN = function(run) {
      micro_state = lift.macro_state(macro_state, parms)
      relaxed_time = time + parms$micro_timestep*parms$micro_relax_steps
      micro_state_relaxed = micro_state_c.stepto(micro_state, parms, control, time = time, timeto = relaxed_time, run = run, record=parms$micro_record)
      next_time = relaxed_time + parms$micro_timestep
      micro_state_next = micro_state_c.stepto(micro_state_relaxed, parms, control, time = relaxed_time, timeto = next_time, run = run, record=parms$micro_record)
      macro_state_relaxed = restrict.micro_state(micro_state_relaxed)
      macro_state_next = restrict.micro_state(micro_state_next)
      return(c(macro_state_relaxed, macro_state_next))
    })
    macro_state_relaxed = rowMeans(vals[1:2,])
    macro_state_next = rowMeans(vals[3:4,])
    macro_state_deriv = (macro_state_next - macro_state)/(2*parms$micro_timestep)
    H = parms$v * macro_state[1] - parms$c * control +
      shadow_state[1] * macro_state_deriv[1] +
      shadow_state[2] * macro_state_deriv[2]
    return(-H)
  }

  opt = nloptr(x0 = control_guess, eval_f = Hamiltonian, lb = parms$control_min, ub = parms$control_max, opts = list(algorithm = "NLOPT_LN_BOBYQA", xtol_rel = 1e-4, xtol_abs=1e-4), macro_state=macro_state, parms=parms, shadow_state=shadow_state, time=time)

  vals = sapply(X = 1:parms$n_sims, FUN = function(run) {
    micro_state = lift.macro_state(macro_state, parms)
    relaxed_time = time + parms$micro_timestep*parms$micro_relax_steps
    micro_state_relaxed = micro_state_c.stepto(micro_state, parms, control = opt$solution, time = time, timeto = relaxed_time, run = run, record=parms$micro_record)
    next_time = relaxed_time + parms$micro_timestep
    micro_state_next = micro_state_c.stepto(micro_state_relaxed, parms, control = opt$solution, time = relaxed_time, timeto = next_time, run = run, record=parms$micro_record)
    macro_state_relaxed = restrict.micro_state(micro_state_relaxed)
    macro_state_next = restrict.micro_state(micro_state_next)
    return(c(macro_state_relaxed, macro_state_next))
  })
  macro_state_relaxed = rowMeans(vals[1:2,])
  macro_state_next = rowMeans(vals[3:4,])
  macro_state_deriv = (macro_state_next - macro_state)/(2*parms$micro_timestep)
  H = parms$v * macro_state[1] - parms$c * opt$solution +
    shadow_state[1] * macro_state_deriv[1] +
    shadow_state[2] * macro_state_deriv[2]
  return(list(control = opt$solution, macro_state_relaxed = macro_state_relaxed, macro_state_next = macro_state_next, macro_state_deriv = macro_state_deriv, hamiltonian = H))
}

alt_shadow_derivs_calc = function(parms=parms, macro_state, macro_deriv, last_deriv_est, macro_second_deriv, shadow_state, control, diff_step) {
  macro_state_alt = macro_state
  macro_state_alt[macro_deriv == 0] = macro_state_alt[macro_deriv == 0] + diff_step

  project_timestep = parms$macro_timestep - parms$micro_timestep*(parms$micro_relax_steps)

  vals = sapply(X = 1:parms$n_sims, FUN = function(run) {
    micro_state = lift.macro_state(macro_state_alt, parms)
    relaxed_time = time + parms$micro_timestep*parms$micro_relax_steps
    micro_state_relaxed = micro_state_c.stepto(micro_state, parms, control, time = time, timeto = relaxed_time, run = run, record=parms$micro_record)
    next_time = relaxed_time + parms$micro_timestep
    micro_state_next = micro_state_c.stepto(micro_state_relaxed, parms, control, time = relaxed_time, timeto = next_time, run = run, record=parms$micro_record)
    macro_state_relaxed = restrict.micro_state(micro_state_relaxed)
    macro_state_next = restrict.micro_state(micro_state_next)
    return(c(macro_state_relaxed, macro_state_next))
  })
  macro_state_alt_relaxed = rowMeans(vals[1:2,])
  macro_state_alt_next = rowMeans(vals[3:4,])
  macro_deriv_alt_forward = (macro_state_alt_next - macro_state)/(2*parms$micro_timestep)


  if(!is.null(last_deriv_est)) {
    macro_deriv_alt = (project_timestep*macro_deriv_alt_forward + (parms$micro_relax_steps*parms$micro_timestep)*last_deriv_est)/parms$macro_timestep
  } else {
    macro_deriv_alt = macro_deriv_alt_forward
  }

  dfdX = (macro_deriv_alt - macro_deriv)/(macro_state_alt - macro_state)


  alt_shadow_derivs = c()
  alt_shadow_derivs[1] = -(parms$v +
                             shadow_state[1] * dfdX[1] +
                             shadow_state[2] * dfdX[1])
  alt_shadow_derivs[2] = -(shadow_state[1] * dfdX[2] +
                             shadow_state[2] * dfdX[2])

  return(alt_shadow_derivs)

}

#' @import polynom
second_deriv_from_3pts = function(x, y) {
  p = poly.calc(x, y)
  if(is.list(p)) {
    dd =2*sapply(p, function(x) `[`(x,3))
  } else {
    dd = 2*p[3]
  }
  dd[is.na(dd)] = 0
  return(dd)
}
