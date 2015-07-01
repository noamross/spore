a =  capture.output(b <- nloptr(x0 = control_guess, eval_f = Hamiltonian, lb = parms$control_min, ub = parms$control_max, opts = list(algorithm = "NLOPT_LN_SBPLX", xtol_rel = 1e-3, print_level = 3), macro_state=macro_state, parms=parms, shadow_state=shadow_state, time=time))
a = paste(a, collapse="\n")
library(stringi)
xvals = as.numeric(stri_extract_all_regex(a, "(?<=x = )[\\d.]+(?=[\\n\\Z]?)")[[1]])
yvals = as.numeric(stri_extract_all_regex(a, "(?<=f\\(x\\) = )[\\d.]+(?=[\\n\\Z]?)")[[1]])
plot(vv
     lines(xvals, yvals)
     points(tail(xvals, 1), tail(yvals,1), col="red", pch=16)

     #' @import stringi magrittr rlist dplyr
     nloptr_trace = function( x0,
                              eval_f,
                              eval_grad_f = NULL,
                              lb = NULL,
                              ub = NULL,
                              eval_g_ineq = NULL,
                              eval_jac_g_ineq = NULL,
                              eval_g_eq = NULL,
                              eval_jac_g_eq = NULL,
                              opts = list(),
                              ...) {

       opts$print_level = 3

       printed_output <- capture.output(out <- nloptr( x0,
                                                       eval_f,
                                                       eval_grad_f,
                                                       lb,
                                                       ub,
                                                       eval_g_ineq,
                                                       eval_jac_g_ineq,
                                                       eval_g_eq,
                                                       eval_jac_g_eq,
                                                       opts,
                                                       ...))

       trace = paste(printed_output, collapse = "\n")
       iterations = trace %>%
         stri_extract_all_regex("(?<=iteration\\:\\s)\\d+(?=\\n)", simplify=TRUE) %>%
         as.vector()

       xvals = trace %>%
         stri_extract_all_regex("(?<=x = \\()[\\s\\-\\d\\,\\.]+", simplify=TRUE) %>%
         stri_trim_both() %>%
         stri_split_regex("[\\s\\,]+", simplify=TRUE)

       fvals = trace %>%
         stri_extract_all_regex("(?<=f\\(x\\) = )[\\s\\-\\d\\,\\.]+", simplify=TRUE) %>%
         stri_trim_both()

       gvals = trace %>%
         stri_extract_all_regex("(?<=g\\(x\\) = )[\\s\\-\\d\\,\\.]+", simplify=TRUE) %>%
         stri_trim_both()

       trace_data = cbind(iterations, xvals, fvals) %>%
         as.data.frame(stringAsFactors=FALSE) %>%
         mutate_each(funs(as.numeric(as.character(.))))

       attr(out, "trace") <- list(iterations, xvals, fvals)
       class(out) <- c(class(out), "traced")
       return(out)
     }


