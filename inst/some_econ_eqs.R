#Returns a single scalar of profit per unit time
profit = function(macro_state, parms, control) {
  macro_state[1]*parms$v - control*parms$c
}

# Returns the derivative of profit with respect to each state variable
profit_derivs = function(macro_state, parms, control) {
  c(parms$v, 0)
  }
#
# Hamiltonian = function(macro_state, parms, control, macro_derivs, shadow_state) {
#   profit(macro_state, parms, control) + sum(shadow_state * macro_derivs)
# }

shadow_derivs = function(macro_state, parms, control, macro_second_derivs, shadow_state) {
  profit_derivs(macro_state, parms, control) + rowSums( (1 / macro_derivs) %o% (shadow_state*macro_second_derivs)
}


