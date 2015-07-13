#' @import nloptr
#' @export
determine_control = function(macro_state, parms, shadow_state, time, control_guess, verbose=FALSE) {

  if(parms$control_min == parms$control_max) {
    opt = parms$control_min
  } else {
  opt = nloptr(x0 = control_guess, eval_f = Hamiltonian, lb = parms$control_min,
               ub = parms$control_max, opts = list(algorithm = "NLOPT_LN_BOBYQA",
                                                   xtol_rel = 1e-3, xtol_abs=1e-3,
                                                   print_level = ifelse(verbose, 3, 0)),
               macro_state=macro_state, parms=parms, shadow_state=shadow_state,
               time=time)$solution
  }

  output = sims_step(macro_state, parms, shadow_state, time, opt)
  return(output)
}


Hamiltonian = function(control, macro_state, parms, shadow_state, time) {
  output = sims_step(macro_state, parms, shadow_state, time, control = control)
  return(-output$hamiltonian)
}

sims_step = function(macro_state, parms, shadow_state, time, control) {
  macro_weights = ceiling(macro_state) - macro_state
  macro_weights[macro_weights == floor(macro_weights)] = 1
  macro_weights_long = apply(expand.grid(lapply(macro_weights, function(z) c(z, 1 -z))), 1, prod)
  macro_state_mod = ifelse(macro_state == floor(macro_state), macro_state + 0.5, macro_state)
  macro_integer_expanded = expand.grid(lapply(macro_state_mod, function(z) c(floor(z), ceiling(z))))
  integer_dims = attr(macro_integer_expanded, "out.attrs")$dim
  macro_integer_list = as.data.frame(t(macro_integer_expanded))

  integer_vals = lapply(macro_integer_list, function(macro_integers) {
    vals = sapply(X = 1:parms$n_sims, FUN = function(run) {
      micro_state = lift_macro_state(macro_integers, parms)
      relaxed_time = time + parms$micro_timestep*parms$micro_relax_steps
      micro_state_relaxed = micro_state_c_stepto(micro_state, parms, control = control,
                                                 time = time, timeto = relaxed_time,
                                                 run = run, record=parms$micro_record)
      next_time = relaxed_time + parms$micro_timestep
      micro_state_next = micro_state_c_stepto(micro_state_relaxed, parms, control = control,
                                              time = relaxed_time, timeto = next_time,
                                              run = run, record=parms$micro_record)
      macro_state_relaxed = restrict.micro_state(micro_state_relaxed)
      macro_state_next = restrict.micro_state(micro_state_next)
      return(c(macro_state_relaxed, macro_state_next))
    })
    return(vals)
  })

  weighted_vals = lapply(1:4, function(i) integer_vals[[i]]*macro_weights_long[i]*4)
  macro_state_relaxed = rowMeans(do.call(cbind, weighted_vals)[1:2,])
  macro_state_next = rowMeans(do.call(cbind, weighted_vals)[3:4,])
  macro_state_deriv = rowMeans(do.call(cbind, weighted_vals)[3:4,] - do.call(cbind, weighted_vals)[1:2,])/parms$micro_timestep
  H = parms$v * macro_state[1] - parms$c * control +
    shadow_state[1] * macro_state_deriv[1] +
    shadow_state[2] * macro_state_deriv[2]

  combos = list(c(1,3), c(2,4), c(1,2), c(3,4))
  diff_weights = list(c(macro_weights[2], 1 - macro_weights[2]),
                      c(macro_weights[2], 1 - macro_weights[2]),
                      c(macro_weights[1], 1 - macro_weights[1]),
                      c(macro_weights[1], 1 - macro_weights[1]))
  integer_derivs = lapply(integer_vals, function(v) {
    rowMeans(v[3:4, ] - v[1:2, ])/parms$micro_timestep
  })

  half_weight_derivs = lapply(1:4, function(i) {
    Reduce(function(v1, v2) {
      (v1*diff_weights[[i]][1] + v2*diff_weights[[i]][2])
    }, integer_derivs[combos[[i]]])
  })

  diff_indices = list(c(1,2), c(3,4))
  macro_dfdx_0 = lapply(diff_indices, function(diff_i) {Reduce(function(v1, v2) v2 - v1, half_weight_derivs[diff_i])})
  macro_dfdx = do.call(rbind, macro_dfdx_0)
  output = list(control = control, macro_state_relaxed = macro_state_relaxed, macro_state_next = macro_state_next, macro_state_deriv = macro_state_deriv, macro_dfdx= macro_dfdx, hamiltonian = H)
  return(output)
}
