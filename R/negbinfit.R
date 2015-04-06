KLD <- function(par, dat) {
  negbin <- dnbinom(x=0:(length(dat)-1), mu=par[1], size=par[2])
  KLD <- sum(dat*pmax(log(dat/negbin), 0))
  return(KLD)
}

NegBinFit <- function(dat, mu.start=NULL, size.start=NULL) {
  dat = dat[1:max(which(dat != 0))]
  dat = dat/sum(dat)
  bins = 0:(length(dat) - 1)
  if (is.null(mu.start)) mu.start = weighted.mean(bins, dat)
  if (is.null(size.start)) {
  	size.start = mu.start^2 / (max(0.01, sum(dat * (bins - mu.start)^2) - mu.start))
  	size.start = ifelse(size.start==0, 1, size.start)
  }
  mu.start = ifelse(mu.start== 0, 0.00001, mu.start)
  return(optim(par=c(mu=mu.start, size=size.start), fn=KLD, dat=dat))
}

simneg = function(n, mu, size) {
  dist = rnbinom(n=n, mu=mu, size=size)
  
}

PoisFit = function(x, lambda_start=NULL) {
	if(is.null(lambda_start)) {
		lambda_start = (x[,1] %*% x[,2]) / sum(x[,2])
	}
	
	result = optim(lambda_start,
								 function(p) {return(-PoisLogLik(x, p))},
			           method = "L-BFGS-B", lower=1e-10)
	
	lambda = result$par
	
	fit = list(lambda = lambda,
						 fitted.values = sum(x[,2]) * dpois(x[,1], lambda),
						 log.likelihood = PoisLogLik(x, lambda))
	return(fit)
}


PoisLogLik = function(x, lambda) {
   if(lambda < 0) {
   	return(-Inf)
   } else {
   	return(x[, 2] %*% dpois(x[, 1], lambda, log = TRUE))
   }
}
