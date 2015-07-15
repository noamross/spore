#' @import nloptr tracer
#' @export
determine_control2 = function(macro_state, parms, shadow_state, time, control_guess, verbose=FALSE, return_trace = FALSE, nloptr_options = NULL) {

  if(is.null(nloptr_options)) {
    nloptr_options = list(
      algorithm = "NLOPT_LN_BOBYQA",
      xtol_rel = 1e-4, xtol_abs=1e-4,
      print_level = ifelse(verbose, 3, 0)
    )}

  if(is.null(nloptr_options$print_level)) {
    nloptr_options$print_level = ifelse(verbose, 3, 0)
  }
  opt = list()
  if(parms$control_min == parms$control_max) {
    opt$solution = parms$control_min
  } else if (!return_trace) {
    opt = nloptr(x0 = control_guess, eval_f = Hamiltonian, lb = parms$control_min,
                 ub = parms$control_max, opts = nloptr_options,
                 macro_state=macro_state, parms=parms, shadow_state=shadow_state,
                 time=time)
  } else {
    opt = nloptr_tr(x0 = control_guess, eval_f = Hamiltonian, lb = parms$control_min,
                    ub = parms$control_max, opts = nloptr_options,
                    macro_state=macro_state, parms=parms, shadow_state=shadow_state,
                    time=time)
  }

  output = sims_step(macro_state, parms, shadow_state, time, opt$solution)
  if(return_trace) {
    output$opt = opt
  }
  return(output)
}


Hamiltonian = function(control, macro_state, parms, shadow_state, time) {
  output = sims_step(macro_state, parms, shadow_state, time, control = control)
  return(-output$hamiltonian)
}

#' @export
sims_step2 = function(macro_state, parms, shadow_state, time, control) {
  macro_weights = ceiling(macro_state) - macro_state
  macro_weights[macro_weights == floor(macro_weights)] = 1
  macro_weights_long = apply(expand.grid(lapply(macro_weights, function(z) c(z, 1 -z))), 1, prod)
  macro_state_mod = ifelse(macro_state == floor(macro_state), macro_state + 0.5, macro_state)
  macro_integer_expanded = expand.grid(lapply(macro_state_mod, function(z) c(floor(z), ceiling(z))))
  integer_dims = attr(macro_integer_expanded, "out.attrs")$dim
  macro_integer_list = as.data.frame(t(macro_integer_expanded))

  integer_vals = lapply(macro_integer_list, function(macro_integers) {
    if(parms$parallel_cores == 1) {
      micro_state = lift_macro_state(macro_integers, parms)
      macro_diff = c(0, 0)
      addtime = 0
      for(i in 1:parms$n_sims) {
        micro_state_next = micro_state_c_step(micro_state, parms, control = control, time = time)
        macro_diff = macro_diff + restrict.micro_state(micro_state_next$micro_state) - macro_integers
        addtime = addtime + micro_state_next$time_next
      }
      macro_diff/addtime
#       vals = sapply(X = 1:parms$n_sims, FUN = function(run) {
#         micro_state_next = micro_state_c_step(micro_state, parms, control = control,
#                                                 time = time)
#         macro_state_rate = (restrict.micro_state(micro_state_next$micro_state) -
#                                restrict.micro_state(micro_state))/((micro_state_next$time_next - time))
#         return(c(restrict.micro_state(micro_state_next$micro_state) - macro_integers,  micro_state_next$time_next - time))
#       })
    } else {
      micro_state = lift_macro_state(macro_integers, parms)
      vals = mclapply(X = 1:parms$n_sims, FUN = function(run) {
        micro_state_next = micro_state_c_step(micro_state, parms, control = control,
                                              time = time)
        macro_state_rate = (restrict.micro_state(micro_state_next$micro_state) -
                              restrict.micro_state(micro_state))/((micro_state_next$time_next - time))
        return(c(restrict.micro_state(micro_state_next$micro_state) - macro_integers,  micro_state_next$time_next - time))
      }, mc.preschedule = TRUE, mc.cores = parms$parallel_cores, mc.silent=TRUE)
      vals = simplify2array(vals)
    }
    return(vals)
  })

  integer_derivs = lapply(integer_vals, function(v) {
    rowMeans(v[1:2,])/mean(v[3,])
  })

  macro_state_deriv = as.vector(simplify2array(integer_derivs) %*% macro_weights_long)

  H = parms$v * macro_state[1] - parms$c * control +
    shadow_state[1] * macro_state_deriv[1] +
    shadow_state[2] * macro_state_deriv[2]

  combos = list(c(1,3), c(2,4), c(1,2), c(3,4))
  diff_weights = list(c(macro_weights[2], 1 - macro_weights[2]),
                      c(macro_weights[2], 1 - macro_weights[2]),
                      c(macro_weights[1], 1 - macro_weights[1]),
                      c(macro_weights[1], 1 - macro_weights[1]))


  half_weight_derivs = lapply(1:4, function(i) {
    Reduce(function(v1, v2) {
      (v1*diff_weights[[i]][1] + v2*diff_weights[[i]][2])
    }, integer_derivs[combos[[i]]])
  })

  diff_indices = list(c(1,2), c(3,4))
  macro_dfdx_0 = lapply(diff_indices, function(diff_i) {Reduce(function(v1, v2) v2 - v1, half_weight_derivs[diff_i])})
  macro_dfdx = do.call(rbind, macro_dfdx_0)
  output = list(control = control, macro_state_deriv = macro_state_deriv, macro_dfdx= macro_dfdx, hamiltonian = H)
  return(output)
}

