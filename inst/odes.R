# ODE problems

library(deSolve)
library(bvpSolve)
library(rootSolve)
library(nloptr)

profout = function(out) {
  integrate(splinefun(out[,"time"], (out[,"N"]*parms$v - (out[,"h"])*parms$cost)*exp(-out[,"time"]*parms$delta)),
            lower = out[1,"time"], upper=tail(out[,"time"], 1))$value
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
    h_calc = uniroot.all(hH_h,c(parms$control_min, parms$control_max), n=10000, t=t, state=state, parms=parms)
    H_calc = H(h_calc, t, state, parms)
    if(length(h_calc) == 1) h_calc = c(h_calc, NA)
    h = max(min(max(h_calc, na.rm=TRUE), control_max), control_min)
    dN = max(r*N*(1 - N/K), 0) - alpha*P - d*N
    if(isTRUE(all.equal(dN, 0))) dN = 0
    dP = lambda*P*N - mu*P - d*P - alpha*P - alpha*(P^2)/N + exp(-h)*N*lambda_ex
    dS1 = -v - S1*(r - d - 2*(r/K)*N) - S2*(lambda*P + alpha*(P^2/N^2) + exp(-h) * lambda_ex) + delta*S1
    dS2 = S1*alpha - S2*(lambda*N - mu - d - alpha - 2*alpha*P/N) + delta*S2
    derivs = c(dN=dN, dP=dP, dS1=dS1, dS2=dS2)
    #    derivs = ifelse(abs(derivs) > .Machine$double.eps ^ 0.5, derivs, 0)
    return(list(derivs, c(derivs, h_calc = h_calc, H_calc = H_calc, h=h)))
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

maxlim = max(as.vector(cbind(out[,c("N", "P","h")], out2[,c("N", "P","h")], out0[,c("N", "P")], out3[,c("N", "P", "h")])))

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
opt0 <- nloptr_tr(x0 = c(0,0), eval_f = profun, lb=c(0, -1e5), ub=c(1e5, 0),
               opts = list(algorithm = "NLOPT_GN_CRS2_LM", maxeval=250, print_level=3), sys=syseq, parms=parms)
opt0

opt0trace = tracer(opt0)
opt = nloptr_tr(x0 = opt0$solution, eval_f = profun, lb=c(-Inf, -Inf), ub=c(Inf, Inf),
             opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval=1000, print_level=3), sys=syseq, parms=parms)

opttrace = tracer(opt)

optshoot0 =  nloptr(x0 = unname(init[3:4]), eval_f = shootfun, lb=c(0, -1e5), ub=c(1e5, 0),
                   opts = list(algorithm = "NLOPT_GN_CRS2_LM", maxeval=500, print_level=3), sys=syseq, parms=parms, target=c(0,0))
optshoot = nloptr(x0 = optshoot0$solution, eval_f = shootfun, lb=c(-Inf, -Inf), ub=c(Inf, Inf),
                 opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval=1000, print_level=3), sys=syseq, parms=parms, target=c(0,0))

# opt2 = nloptr(x0 = unname(c(1933, -3261)), eval_f = profun, lb=c(-Inf, -Inf), ub=c(Inf, 0),
#               opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval=10000), sys=syseq2, parms=parms)
# opt2


#  Scenario based on SIR - adding more hosts
#  Add the discount?
#  Write a crude outline of the paper. (Author list?)
#
