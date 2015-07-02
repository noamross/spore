# ODE problems

library(deSolve)
library(bvpSolve)
library(rootSolve)
library(nloptr)

profout = function(out) {
  intval = integrate(splinefun(out[,"time"], (out[,"N"]*parms$v - (out[,"h"])*parms$cost)*exp(-out[,"time"]*parms$delta)),
            lower = out[1,"time"], upper=tail(out[,"time"], 1))
  return(intval$value)
}

hH_h = function(h, t, state, parms) {
  with(as.list(c(state, parms)), {
    h*( - cost - S2 * exp(-h) * N * lambda_ex)
  })
}

H = function(h, t, state, parms) {
  with(as.list(c(state, parms)), {
    val = v*N - cost*h + S1*(r*N*(1 - N/K) - alpha*P - d*N) + S2*(lambda*P*N - mu*P - d*P - alpha*P - alpha*(P^2)/N + exp(-h)*N*lambda_ex)
    return(val)
  })
}

# control = function(h, t, state, parms, rule) {
#  h = uniroot.all(hH_h,c(parms$control_min, parms$control_max), n=10000, t=t, state=state, parms=parms)
#  switch(rule,
#    max =
#  )
#   }

syseq = function(t, state, parms) {
  with(as.list(c(state, parms)), {
    h_calc = suppressWarnings(log(- S2 * N * lambda_ex/cost))
    if(is.nan(h_calc)) h = 0 else h = h_calc
    h = min(max(h, control_min), control_max)
    H_calc = H(h, t, state, parms)
    dN = max(r*N*(1 - N/K), 0) - alpha*P - d*N
    if(isTRUE(all.equal(dN, 0))) dN = 0
    dP = P*(lambda*N - mu - d - alpha - alpha*(P)/N) + exp(-h)*N*lambda_ex
    dS1 = -v - S1*(r - d - 2*(r/K)*N) - S2*(lambda*P + alpha*(P^2/N^2) + exp(-h) * lambda_ex) + delta*S1
    dS2 = S1*alpha - S2*(lambda*N - mu - d - alpha - 2*alpha*P/N) + delta*S2
    derivs = c(dN=dN, dP=dP, dS1=dS1, dS2=dS2)
    return(list(derivs, c(derivs, h_calc=h_calc, H_calc = H_calc, h=h)))
  })
}

syseq2 = function(t, state, parms) {
  with(as.list(c(state, parms)), {
    h_calc = c(parms$control_min, parms$control_max, uniroot.all(hH_h,c(parms$control_min, parms$control_max), n=10000, t=t, state=state, parms=parms))
    if(length(h_calc) == 3) h_calc = c(h_calc, NA)
    H_calc = H(h_calc, t, state, parms)
    if(length(H_calc) == 3) H_calc = c(h_calc, NA)
    H_sel = max(H_calc, na.rm=TRUE)
    h = h_calc[which.max(H_calc)]
 #   browser()
    dN = max(r*N*(1 - N/K), 0) - alpha*P - d*N
    if(isTRUE(all.equal(dN, 0))) dN = 0
    dP = lambda*P*N - mu*P - d*P - alpha*P - alpha*(P^2)/N + exp(-h)*N*lambda_ex
    dS1 = -v - S1*(r - d - 2*(r/K)*N) - S2*(lambda*P + alpha*(P^2/N^2) + exp(-h) * lambda_ex) + delta*S1
    dS2 = S1*alpha - S2*(lambda*N - mu - d - alpha - 2*alpha*P/N) + delta*S2
    derivs = c(dN=dN, dP=dP, dS1=dS1, dS2=dS2)
    #    derivs = ifelse(abs(derivs) > .Machine$double.eps ^ 0.5, derivs, 0)
    return(list(derivs, c(derivs, h_calc = h_calc, H_calc = H_calc, H_sel=H_sel, h=h)))
  })
}


syseq3 = function(t, state, parms) {
  with(as.list(c(state, parms)), {
    h_calc = suppressWarnings(log(- S2 * N * lambda_ex/cost))
    if(is.nan(h_calc)) {
      h = 0
      H_calc = c(NA, H(h, t, state, parms))
    } else {
      H_calc = H(c(0, h_calc), t, state, parms)
      h = c(0, h_calc)[which.max(H_calc)]
    }
    dN = max(r*N*(1 - N/K), 0)  - alpha*P - d*N
    if(isTRUE(all.equal(dN, 0))) dN = 0
    dP = lambda*P*N - mu*P - d*P - alpha*P - alpha*(P^2)/N + exp(-h)*N*lambda_ex
    dS1 = -v - S1*(r - d - 2*(r/K)*N) - S2*(lambda*P + alpha*(P^2/N^2) + exp(-h)*lambda_ex) + delta*S1
    dS2 = S1*alpha - S2*(lambda*N - mu - d - alpha - 2*alpha*P/N) + delta*S2
    derivs = c(dN=dN, dP=dP, dS1=dS1, dS2=dS2)
    #    derivs = ifelse(abs(derivs) > .Machine$double.eps ^ 0.5, derivs, 0)
    return(list(derivs, c(derivs, h_calc = h_calc, H_calc= H_calc, h=h)))
  })
}


parms = list(
  lambda = 0.001,
  lambda_ex = 0.2,
  alpha = 0.1,
  mu = 0.01,
  r = 0.5,
  d = 0.01,
  K = 100,
  control_min = 0,
  control_max = 1000,
  v = 50,
  cost = 200,
  delta = 0,
  progress = TRUE
  #micro_record = file("micro.txt", open="w+")
  #  macro_record = file("macro.txt", open="w")
)


N0 = 100
P0 = 0
#S10 = with(parms, {-v/(r - d - 2*(r/K)*N0)})
S10 = 97.323969
#S20 = with(parms, {(S10*alpha)/(lambda*N0 - mu - d - alpha - 2*alpha*P0/N0)})
S20 =  -306.152002
init = c(N = N0, P = P0, S1 = S10, S2 = S20)
times = seq(0,40, by=1)

parms0 = parms
parms0$control_max = 0

out0 = lsoda(init, times, syseq, parms0)
out = lsoda(init, times, syseq, parms)
out2 = lsoda(init, times, syseq2, parms)
out3 = lsoda(init, times, syseq3, parms)

#maxlim = max(as.vector(cbind(out[,c("N", "P","h")], out2[,c("N", "P","h")], out0[,c("N", "P")], out3[,c("N", "P", "h")])))
maxlim = max( out2[,c("N", "P","h")])
par(mfrow = c(1,4))
matplot(out0[,1], out0[,c("N", "P", "h")], type="l", ylim=c(0,maxlim), lty=1, lwd=1, main=round(profout(out0)))
matplot(out[,1], out[,c("N", "P", "h")], type="l", ylim=c(0,maxlim), lty=1, lwd=3, main=round(profout(out)))
matplot(out2[,1], out2[,c("N", "P", "h")], type="l", ylim=c(0,maxlim), lty=1, lwd=1, main=round(profout(out2)))
matplot(out3[,1], out3[,c("N", "P", "h")], type="l", ylim=c(0,maxlim), lty=1, lwd=1, main=round(profout(out3)))

# maxlim2 = max(as.vector(cbind(out[,c("S1", "S2", "h")], out2[,c("S1", "S2", "h")], out0[,c("S1", "S2", "h")])))
# minlim2 = min(as.vector(cbind(out[,c("S1", "S2", "h")], out2[,c("S1", "S2", "h")], out0[,c("S1", "S2", "h")])))
# matplot(out0[,1], out0[,c("S1", "S2", "h")], type="l", ylim=c(minlim2,maxlim2), lty=1, lwd=1, main=round(profout(out0)))
# matplot(out[,1], out[,c("S1", "S2", "h")], type="l", ylim=c(minlim2,maxlim2), lty=1, lwd=3, main=round(profout(out)))
# matplot(out2[,1], out2[,c("S1", "S2", "h")], type="l", ylim=c(minlim2,maxlim2), lty=1, lwd=1, main=round(profout(out2)))


shootfun = function(ivals, sys, target, parms) {
  init = c(N=N0, P=P0, S1=ivals[1], S2=ivals[2])
  out = lsoda(init, times, sys, parms)
  obj = sum((tail(out, 1)[,c("S1", "S2")] - target)^2)
  if(any(is.na(obj))) obj = Inf
  return(obj)
}

shootfun2 = function(ivals, sys, target, parms) {
  init = c(N=N0, P=P0, S1=ivals[1], S2=ivals[2])
  out = lsoda(init, times, sys, parms)
  obj = sum((tail(out, 1)[,c("N", "P")] - target)^2)
  if(any(is.na(obj))) obj = Inf
  return(obj)
}

optfun = function(ivals, sys, parms) {
  init = c(N=N0, P=P0, S1=ivals[1], S2=ivals[2])
  out = lsoda(init, times, sys, parms)
  # obj = sum(out[,"h"])
  obj = sum((tail(out, 1)[4:5] - c(0, 0))^2)
  if(any(is.na(obj))) obj = Inf
  return(obj)

}

profun = function(ivals, sys, parms) {
  init = c(N=N0, P=P0, S1=ivals[1], S2=ivals[2])
  out = lsoda(init, times, sys, parms)
  -profout(out)
}


library(tracer)
# #
opt0 <- nloptr_tr(x0 = c(S10,S20), eval_f = profun, lb=c(-1e5, -1e5), ub=c(1e5, 1e5),
               opts = list(algorithm = "NLOPT_GN_CRS2_LM", maxeval=250, print_level=3), sys=syseq2, parms=parms)
opt0

opt0trace = tracer(opt0)
opt = nloptr_tr(x0 = opt0$solution, eval_f = profun, lb=c(-Inf, -Inf), ub=c(Inf, Inf),
             opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval=1000, print_level=3), sys=syseq2, parms=parms)

opttrace = tracer(opt)

out_opt = lsoda(c(N = N0, P = P0, S1 = opt$solution[1], S2 = opt$solution[2]), times, syseq2, parms)

matplot(out_opt[,1], out_opt[,c("N", "P", "h")], type="l", lty=1, lwd=1, main=round(profout(out_opt)))

targets = unname(tail(out_opt, 1)[, c("S1", "S2")])

optshoot0 =  nloptr_tr(x0 = c(0,0), eval_f = shootfun, lb=c(-1e5, -1e5), ub=c(1e5, 1e5),
                   opts = list(algorithm = "NLOPT_GN_CRS2_LM", maxeval=250, print_level=3), sys=syseq2, parms=parms, target=targets)
optshoot0_tr = tracer(optshoot0)
optshoot = nloptr_tr(x0 = optshoot0$solution, eval_f = shootfun, lb=c(-Inf, -Inf), ub=c(Inf, Inf),
                 opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval=1000, print_level=3), sys=syseq2, parms=parms, target=targets)
optshoot2 = nloptr_tr(x0 = c(optshoot0$solution[1], -465), eval_f = shootfun, lb=c(-Inf, -Inf), ub=c(Inf, Inf),
                 opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval=1000, print_level=3), sys=syseq2, parms=parms, target=targets)

out_shoot = lsoda(c(N = N0, P = P0, S1 = optshoot$solution[1], S2 = optshoot$solution[2]), times, syseq2, parms)
profout(out_shoot)
matplot(out_shoot[,1], out_shoot[,c("N", "P", "h")], type="l", lty=1, lwd=1, main=round(profout(out_shoot)))

targets2 = unname(tail(out_opt, 1)[, c("N", "P")])

optshoot_state0 = nloptr_tr(x0 = c(0,0), eval_f = shootfun2, lb=c(-1e5, -1e5), ub=c(1e5, 1e5),
                   opts = list(algorithm = "NLOPT_GN_CRS2_LM", maxeval=250, print_level=3), sys=syseq2, parms=parms, target=targets2)

optshoot_state = nloptr_tr(x0 = c(optshoot_state0$solution[1], -500), eval_f = shootfun2, lb=c(-Inf, -Inf), ub=c(Inf, Inf),
                 opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval=1000, print_level=3), sys=syseq2, parms=parms, target=targets2)

optshoot_st_tr = tracer(optshoot_state0)

out_shoot_st = lsoda(c(N = N0, P = P0, S1 = optshoot_state$solution[1], S2 = optshoot_state$solution[2]), times, syseq2, parms)
matplot(out_shoot_st[,1], out_shoot_st[,c("N", "P", "h")], type="l", lty=1, lwd=1, main=round(profout(out_shoot)))

plot3d(optshoot_st_tr$xval1, optshoot_st_tr$xval2, log(optshoot_st_tr$fval))

# opt2 = nloptr(x0 = unname(c(1933, -3261)), eval_f = profun, lb=c(-Inf, -Inf), ub=c(Inf, 0),
#               opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval=10000), sys=syseq2, parms=parms)
# opt2


#  Scenario based on SIR - adding more hosts
#  Add the discount?
#  Write a crude outline of the paper. (Author list?)
#
