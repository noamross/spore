#devtools::install_github('noamross/cmp')
library(cmp)
library(rethinking)
library(rstan)
dat = data.frame(counts = rcmp(100, 2, 5), counts2 = rcmp(100, 2, 0.5))

model1_pois = alist(
  counts ~ dpois(lambda) ,
  log(lambda) <- a,
  a ~ dnorm(0,10))

fit1_pois = map(model1_pois, data=dat)

model1_cmp = alist(
  counts ~ dcmp(lambda, nu),
  log(lambda) <- a,
  log(nu) <- b,
  a ~ dnorm(0,10),
  b ~ dnorm(0,10))

model1_nb = alist(
  counts ~ dnbinom(mu=lambda, size=nu),
  log(lambda) <- a,
  log(nu) <- b,
  a ~ dnorm(0,10),
  b ~ dnorm(0,10))


fit1_cmp = map(model1_cmp, data=dat)

model2_pois = alist(
  counts2 ~ dpois(lambda),
  log(lambda) <- a,
  a ~ dnorm(0,10))

model2_cmp = alist(
  counts2 ~ dcmp(lambda, nu),
  log(lambda) <- a,
  log(nu) <- b,
  a ~ dnorm(0,10),
  b ~ dnorm(0,10))

model2_nb = alist(
  counts2 ~ dnbinom(mu=lambda, size=nu),
  log(lambda) <- a,
  log(nu) <- b,
  a ~ dnorm(0,10),
  b ~ dnorm(0,10))


fit1_pois = map(model1_pois, data=dat)
fit1_cmp = map(model1_cmp, data=dat)
fit1_nb = map(model1_nb, data=dat)
fit2_pois = map(model2_pois, data=dat)
fit2_cmp = map(model2_cmp, data=dat)
fit2_nb = map(model2_nb, data=dat)


compare(fit1_cmp, fit1_pois, fit1_nb, WAIC=FALSE)
compare(fit2_cmp, fit2_pois, fit2_nb, WAIC=FALSE)

compare(fit1_cmp, fit1_pois, fit1_nb, WAIC=TRUE)
compare(fit2_cmp, fit2_pois, fit2_nb, WAIC=TRUE)

