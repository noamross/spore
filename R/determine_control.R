save_guess = new.env()

#' @export
opt_derivs = function(t, y, parms) {
  parms = relist(parms)
  if(is.null(save_guess$guess)) save_guess$guess = 5
  derivs = determine_control(y[1:2], parms, y[3:4], t,save_guess$guess)
  save_guess$guess = derivs$control

  dS1 = -parms$v - y[3] * derivs$macro_state_jacobian[1, 1] - y[4] * derivs$macro_state_jacobian[1, 2]
  dS2 = -y[3] * derivs$macro_state_jacobian[2, 1] - y[4] * derivs$macro_state_jacobian[2, 2]


  output = c(dN = derivs$macro_state_deriv[1],
             dP = derivs$macro_state_deriv[2],
             dS1 = dS1, dS2 = dS2)
  cat(t, y, output, derivs$control, derivs$H, "\n", sep = "\t")
  cat(t, y, output, derivs$macro_state_jacobian, derivs$control, derivs$H, "\n", file = "out.txt", sep = ",", append = TRUE)
  return(list(output))
}

#' @import nloptr tracer
#' @export
determine_control = function(macro_state, parms, shadow_state, time, control_guess, verbose=FALSE, return_trace = FALSE, nloptr_options = NULL) {

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
    opt = nloptr(x0 = control_guess, eval_f = objective, lb = parms$control_min,
                 ub = parms$control_max, opts = nloptr_options,
                 macro_state=macro_state, parms=parms, shadow_state=shadow_state,
                 time=time)
  } else {
    opt = nloptr_tr(x0 = control_guess, eval_f = objective, lb = parms$control_min,
                    ub = parms$control_max, opts = nloptr_options,
                    macro_state=macro_state, parms=parms, shadow_state=shadow_state,
                    time=time)
  }

  result = macro_state_c_deriv_aves(macro_state, parms, time, opt$solution)
  result$H = Hamiltonian(macro_state, result$macro_state_deriv, shadow_state, parms, time, opt$solution)
  result$control = opt$solution
  if(return_trace) {
    result$opt = opt
  }
  return(result)
}


Hamiltonian = function(macro_state, macro_state_deriv, shadow_state, parms, time, control) {
  return(parms$v * macro_state[1] - parms$c * control +
    shadow_state[1] * macro_state_deriv[1] +
    shadow_state[2] * macro_state_deriv[2])
}

objective = function(control, macro_state, shadow_state, parms, time) {
  derivs = macro_state_c_deriv_aves(macro_state = macro_state, parms = parms, time = time, control = control)
  output = -Hamiltonian(macro_state = macro_state, macro_state_deriv = derivs$macro_state_deriv, shadow_state = shadow_state, parms = parms, time = time, control = control)
  return(output)
}

#' @export
macro_state_c_deriv_aves = function(macro_state, parms, time, control) {
  macro_weights = ceiling(macro_state) - macro_state
  macro_weights[macro_weights == floor(macro_weights)] = 1
  macro_weights_long = apply(expand.grid(lapply(macro_weights, function(z) c(z, 1 -z))), 1, prod)
  macro_state_mod = ifelse(macro_state == floor(macro_state), macro_state + 0.5, macro_state)
  macro_integer_expanded = expand.grid(lapply(macro_state_mod, function(z) c(floor(z), ceiling(z))))
  integer_dims = attr(macro_integer_expanded, "out.attrs")$dim
  macro_integer_list = as.list(as.data.frame(t(macro_integer_expanded)))
  RNGkind("L'Ecuyer-CMRG")
  integer_derivs = mclapply(macro_integer_list, mc.preschedule = FALSE, mc.cores = parms$parallel_cores, mc.set.seed = TRUE,
                            FUN = function(macro_integers) {
                              a = macro_state_c_step_aves(macro_integers, parms, control, time)
                              return(a)
                            })

  weighted_ave_derivs = mapply(`*`, integer_derivs, macro_weights_long, SIMPLIFY = FALSE)
  ave_deriv = Reduce(`+`, weighted_ave_derivs, c(0,0))

  #return(ave_deriv)

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
  output = list(macro_state_deriv = ave_deriv, macro_state_jacobian= macro_dfdx)
  return(output)
}

