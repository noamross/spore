profit = function(macro_state, parms, control) {
  macro_state[1]*parms$v - control*parms$c
}

profit_derivs = function(macro_state, parms, control) {
  c(parms$v, 0)
  }

Hamiltonian = profit(macro_state, parms, control) + shadow_state %*% macro_derivs

shadow_derivs = profit_derivs(macro_state, parms, control) +
  rowSums( (1 / macro_derivs) %o% (shadow_state*macro_second_derivs) )
