a =  capture.output(b <- nloptr(x0 = control_guess, eval_f = Hamiltonian, lb = parms$control_min, ub = parms$control_max, opts = list(algorithm = "NLOPT_LN_SBPLX", xtol_rel = 1e-3, print_level = 3), macro_state=macro_state, parms=parms, shadow_state=shadow_state, time=time))
a = paste(a, collapse="\n")
library(stringi)
xvals = as.numeric(stri_extract_all_regex(a, "(?<=x = )[\\d.]+(?=[\\n\\Z]?)")[[1]])
yvals = as.numeric(stri_extract_all_regex(a, "(?<=f\\(x\\) = )[\\d.]+(?=[\\n\\Z]?)")[[1]])
plot(vv
lines(xvals, yvals)
points(tail(xvals, 1), tail(yvals,1), col="red", pch=16)
