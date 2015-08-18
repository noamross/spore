## ----ODEsoln------------------------------------------------------------


h_fun = function(x) {
  ifelse(x >= 1, log(x), 0)
}

syseq2 = function(t, state, parms) {
  with(as.list(c(state, parms)), {
    if(parms$nocontrol) {
      h = 0
    } else {
      h = h_fun(-S2 * N * lambda_ex / c)
    }
    dN = max(r*N*(1 - N/K), 0) - alpha*P - d*N
    if(isTRUE(all.equal(dN, 0))) dN = 0
    dP = lambda*P*N - mu*P - d*P - alpha*P - alpha*(P^2)/N + exp(-h)*N*lambda_ex
    ddNdN= r - d - 2*(r/K)*N
    ddPdN=lambda*P + alpha*(P^2/N^2) + exp(-h) * lambda_ex
    ddNdP=-alpha
    ddPdP=lambda*N - mu - d - alpha - 2*alpha*P/N
    dS1 = -v - S1*ddNdN - S2*ddPdN + delta*S1
    dS2 = -S1*ddNdP - S2*ddPdP + delta*S2
    derivs = c(dN=dN, dP=dP, dS1=dS1, dS2=dS2)
    return(list(derivs, c(derivs, ddNdN=ddNdN, ddPdN=ddPdN, ddNdP=ddNdP, ddPdP=ddPdP, h=h)))
  })
}

profout = function(out) {
  intval = integrate(splinefun(out[,"time"], (out[,"N"]*parms$v - (out[,"h"])*parms$c)*exp(-out[,"time"]*parms$delta)),
                     lower = out[1,"time"], upper=tail(out[,"time"], 1))
  return(intval$value)
}

profun_base = function(ivals, sys, parms, times) {
  init = c(N=N0, P=P0, S1=ivals[1], S2=ivals[2])
  out = ode(init, times, sys, parms, method="euler")
  -profout(out)
}

shadshoot_base = function(ivals, sys, parms, times, target, inits) {
  init = c(N=inits[1], P=inits[2], S1=ivals[1], S2=ivals[2])
  out = ode(init, times, sys, parms, method="euler")
  # if(any(out[, "S2"] > 0)) return(Inf)
  obj = sum((tail(out, 1)[,c("S1", "S2")] - target)^2)
  if(any(is.na(obj))) obj = Inf
  return(obj)
}

profun = function(ivals, sys, parms, times) {
  memoizedCall(profun_base, ivals=ivals, sys=sys, parms=parms, times=times)
}

shadshoot = function(ivals, sys, parms, times, target, inits) {
  memoizedCall(shadshoot_base, ivals=ivals, sys=sys, parms=parms, times = times, target=target, inits=inits)
}

N0 = 299.9999
P0 = 5.00001
#guess
parms$nocontrol = FALSE

S10 =  19.96184522561706131683
S20 =   -10.62974659561256984830
#opt0_s <- nloptr(x0 = c(S10,S20), eval_f = shadshoot, lb=c(0, -1e5), ub=c(1e5, 0),
#                    opts = list(algorithm = "NLOPT_GN_CRS2_LM", maxeval=10000, xtol_abs = 1e-3 , ftol_rel = -1, xtol_rel = -1, stopval = sqrt(.Machine$double.eps), print_level = ifelse(interactive(), 3, 0), ranseed = 1), sys=syseq2, parms=parms, times=times, target = c(0,0), inits=c(N0, P0))
#
#next_x0 = opt0_s$solution
opt_s <- nloptr(x0 = c(S10,S20), eval_f = shadshoot, lb=c(0, -500), ub=c(500, 0),
                   opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval=10000, xtol_abs = -1 , xtol_rel = -1, ftol_rel = -1, ftol_abs = -1, stopval = sqrt(.Machine$double.eps), print_level = ifelse(interactive(), 3, 0), ranseed = 1), sys=syseq2, parms=parms, times=times, target = c(0,0), inits=c(N0, P0))

out_opt = ode(c(N = N0, P = P0, S1 = opt_s$solution[1], S2 = opt_s $solution[2]), times, syseq2, parms, method="euler")
out_opt = as.data.frame(out_opt)

## ----fig2------------------------------------------------------------

paths = ggplot() +
  geom_line(data = out_opt, mapping = aes(x = time, y = N), col="blue", lwd=1, lty=1) +
  geom_line(data = out_opt, mapping = aes(x = time, y = P), col="red", lwd=1, lty=1) +
  annotate("text", x=c(5,5), y=c(100, 400), label = c("pathogen", "host"),
           color = c("red", "blue"), size=5) +
  labs(list(Title = "", x = "", y = "Population")) + theme_nr

shadpaths = ggplot() +
  geom_line(data = out_opt, mapping = aes(x = time, y = S1), col="darkblue", lwd=1, lty=1) +
  geom_line(data = out_opt, mapping = aes(x = time, y = S2), col="darkred", lwd=1, lty=1) +
  annotate("text", x=c(6,2), y=c(-5, 10), label = c("pathogen shadow value", "host shadow value"),
           color = c("darkred", "darkblue"), size=5) +
  geom_hline(yintercept=0, lty=2, lwd=0.5, col="black") +
  labs(list(Title = "", x = "", y = "Shadow Value")) + theme_nr

controlpath = ggplot() +
  geom_line(data = out_opt, mapping = aes(x = time, y = h), col="black", lwd=1, lty=1) +
  labs(list(Title = "", x = "Time", y = "Control")) + theme_nr

FIG2 = cowplot::plot_grid(plot_grid(paths, shadpaths, ncol=2, labels = c('A', 'B')),
                   controlpath, labels = c("", "C"), ncol = 1)

## ----fig3-----------------------------------------------------------------------

efcontrol = read_csv('out_coarse.txt')
names(efcontrol) = c("time", "N", "P", "S1", "S2", "dN", "dP", "dS1", "dS2", "ddNdN", "ddPdN", "ddNdP", "ddPdP", "h", "H")

paths = ggplot() +
  geom_line(data = efcontrol, mapping = aes(x = time, y = N), col="blue", lwd=1, lty=1) +
  geom_line(data = efcontrol, mapping = aes(x = time, y = P), col="red", lwd=1, lty=1) +
  annotate("text", x=c(5,5), y=c(100, 400), label = c("pathogen", "host"),
           color = c("red", "blue"), size=5) +
  labs(list(Title = "", x = "", y = "Population")) + theme_nr

shadpaths = ggplot() +
  geom_line(data = efcontrol, mapping = aes(x = time, y = S1), col="darkblue", lwd=1, lty=1) +
  geom_line(data = efcontrol, mapping = aes(x = time, y = S2), col="darkred", lwd=1, lty=1) +
  annotate("text", x=c(6,2), y=c(-8, 10), label = c("pathogen shadow value", "host shadow value"),
           color = c("darkred", "darkblue"), size=5) +
  geom_hline(yintercept=0, lty=2, lwd=0.5, col="black") +
  labs(list(Title = "", x = "", y = "Shadow Value")) + theme_nr

controlpath = ggplot() +
  geom_line(data = efcontrol, mapping = aes(x = time, y = h), col="black", lwd=1, lty=1) +
  labs(list(Title = "", x = "Time", y = "Control")) + theme_nr

FIG3 = cowplot::plot_grid(plot_grid(paths, shadpaths, ncol=2, labels = c('A', 'B')),
                   controlpath, labels = c("", "C"), ncol = 1)

## ----fig4-----------------------------------------------------------------------

profits = data.frame(method = factor(c("No Control - ODE", "No Control - EF", "Optimal Control - ODE", "Optimal Control - EF"), levels = c("No Control - ODE", "No Control - EF", "Optimal Control - ODE", "Optimal Control - EF")),
                     profit = c(
                       integrate(approxfun(x = out0$time, y = out0$N*parms$v - out0$h*parms$c), 0, 10)$value,
                       integrate(approxfun(x = as.data.frame(ef_sims)$time, y = as.data.frame(ef_sims)$N*parms$v), 0, 10)$value,
                       integrate(approxfun(x = out_opt$time, y = out_opt$N*parms$v - out_opt$h*parms$c), 0, 10)$value,
                       integrate(approxfun(x = efcontrol$time, y = efcontrol$N*parms$v - efcontrol$h*parms$c), 0, 10)$value
                     ))

FIG4 = ggplot(profits, aes(x=method, y = profit)) + geom_bar(stat="identity") + theme_nr +
  scale_y_continuous(expand=c(0,0)) +
  xlab("") + ylab("Net Profit")


# ## ---ESAplots2-------------------------------------------
#
# ggplot() +
#  	geom_line(data = out_opt, mapping = aes(x = times, y = N), col="steelblue", lwd=3, lty=1) +
#   geom_line(data = out_opt, mapping = aes(x = times, y = P), col="tomato", lwd=3, lty=1) +
#   scale_x_continuous(limits = c(0,10), breaks=c(0, 2.5, 5, 7.5, 10)) +
#   labs(list(Title = "", x = "", y = "")) + theme_nr +
#   theme(axis.text = element_text(size=20))
#
# ggplot() +
#   geom_line(data = out_opt, mapping = aes(x = time, y = h), col="darkgrey", lwd=3, lty=1) +
#   labs(list(Title = "", x = "", y = "")) + theme_nr +
#   theme(axis.text = element_text(size=20))
#
# efcontrol = read_csv('out_coarse.txt')
# names(efcontrol) = c("time", "N", "P", "S1", "S2", "dN", "dP", "dS1", "dS2", "ddNdN", "ddPdN", "ddNdP", "ddPdP", "h", "H")
#
# ggplot() +
#  	geom_line(data = out_opt, mapping = aes(x = times, y = N), col="steelblue", lwd=3, lty=1) +
#   geom_line(data = out_opt, mapping = aes(x = times, y = P), col="tomato", lwd=3, lty=1) +
#   geom_line(data = efcontrol, mapping = aes(x = times, y = N), col="black", lwd=3, lty=3) +
#   geom_line(data = efcontrol, mapping = aes(x = times, y = P), col="black", lwd=3, lty=3) +
#   scale_x_continuous(limits = c(0,10), breaks=c(0, 2.5, 5, 7.5, 10)) +
#   labs(list(Title = "", x = "", y = "")) + theme_nr +
#   theme(axis.text = element_text(size=20))
#
# ggplot() +
#   geom_line(data = out_opt, mapping = aes(x = time, y = h), col="darkgrey", lwd=3, lty=1) +
#   geom_line(data = efcontrol, mapping = aes(x = times, y = h), col="black", lwd=3, lty=3) +
#   labs(list(Title = "", x = "", y = "")) + theme_nr +
#   theme(axis.text = element_text(size=20))
#
#
# ggplot() +
#   geom_line(data = out0, mapping = aes(x = time, y = N), col="steelblue", lwd=2, lty=1) +
#   geom_line(data = out0, mapping = aes(x = time, y = P), col="tomato", lwd=2, lty=1) +
#  	geom_line(data = as.data.frame(ef_sims), mapping = aes(x = times, y = N), col="black", lwd=3, lty=3) +
#   geom_line(data = as.data.frame(ef_sims), mapping = aes(x = times, y = P), col="black", lwd=3, lty=3) +
#   scale_x_continuous(limits = c(0,10), breaks=c(0, 2.5, 5, 7.5, 10)) +
#   labs(list(Title = "", x = "", y = "")) + theme_nr +
#   theme(axis.text = element_text(size=20))
#
# ggplot() +
#   geom_line(data = out0, mapping = aes(x = time, y = N), col="steelblue", lwd=2, lty=1) +
#   geom_line(data = out0, mapping = aes(x = time, y = P), col="tomato", lwd=2, lty=1) +
#  	geom_line(data = as.data.frame(ef_sims), mapping = aes(x = times, y = N), col="black", lwd=3, lty=3) +
#   geom_line(data = as.data.frame(ef_sims), mapping = aes(x = times, y = P), col="black", lwd=3, lty=3) +
#   geom_line(data = macro_output_mean, mapping = aes(x = time, y = mean_N), col="darkblue", lwd = 2, lty=2) +
# 	geom_line(data = macro_output_mean, mapping = aes(x = time, y = mean_P), col="darkred", lwd = 2, lty=2) +
#   scale_x_continuous(limits = c(0,10), breaks=c(0, 2.5, 5, 7.5, 10)) +
#   labs(list(Title = "", x = "", y = "")) + theme_nr +
#   theme(axis.text = element_text(size=20))
#
# ggplot() +
# 	geom_line(data = as.data.frame(ef_sims), mapping = aes(x = times, y = N), col="darkblue", lwd=2, lty=6) +
#   geom_line(data = as.data.frame(ef_sims), mapping = aes(x = times, y = P), col="darkred", lwd=2, lty=6) +
#   scale_x_continuous(limits = c(0,10), breaks=c(0, 2.5, 5, 7.5, 10)) +
#   labs(list(Title = "", x = "", y = "")) + theme_nr +
#   theme(axis.text = element_text(size=20))
#
#
