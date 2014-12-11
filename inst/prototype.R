parms = list(
	max_i = 100,
	lambda = 1,
	lambda_ex = 1,
	alpha = 0.5,
	mu = 0.1,
	r = 0.2,
	d = 0.05,
	K = 100,
	init_pop = 100,
	t_max = 100
)

t = 0
pop = rep(0, parms$max_i + 1)
is = 0:(parms$max_i)
pop[1] = parms$init_pop

fname = "out.txt"
cat(c(t, pop, "\n"), file=fname, append=FALSE)

while (t < parms$t_max) {
	infections = sum(is*pop)
  N = sum(pop)
	if (N == 0) break()
	probs = c(
		P_inf = parms$lambda * infections + parms$lambda_ex*N,
		P_rec = parms$mu * infections,
		P_die = parms$alpha * infections + parms$d*N,
		P_birth = parms$r * N * (1 - N/parms$K)
	)
	t_next = t + rexp(1, probs)
  event = sample(1:4, 1, prob=probs)	
  if (event == 4) {
  	pop[1] = pop[1] + 1
  } else if (event == 3) {
		i = sample(is, 1, prob=parms$alpha * is * pop + parms$d*pop)
		pop[i+1] = pop[i+1] - 1
  } else if (event == 2) {
  	i = sample(is, 1, prob = parms$mu * is * pop)
  	pop[i+1] = pop[i+1] - 1
  	pop[i] = pop[i] + 1
  } else if (event == 1) {
  	i = sample(is, 1, prob=pop)
  	pop[i + 1] = pop[i + 1] - 1
  	pop[i + 2] = pop[i + 2] + 1
  }
	#browser()
	cat(c(t_next, pop, "\n"), file=fname, append=TRUE)
	t = t_next
}

output = read.table(fname, sep = " ")
colnames(output) = c("time", 0:100)
output[, parms$max_i + 3] = NULL
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggvis)
df = output %>%
	gather(infections, population, -time, convert=TRUE) %>%
	arrange(time, infections)

# df %>%
#  ggvis(~infections, ~population) %>%
#  filter(time == max(time >= eval(input_slider(0, 0.05)))) %>%
#  layer_bars()

dfv = df %>%
	group_by(time) %>%
	summarize(v = max(infections[population > 0]))
max(dfv$v, na.rm=TRUE)

# manipulate(
# 	ggplot(subset(df, time == max(df$time[df$time <= TIME])),
# 				 aes(x=infections, y=population)) + 
# 		geom_bar(stat="identity") + 
# 		ylim(c(0, max(df$population))) +
# 		xlim(c(-0.5, max(dfv$v, na.rm=TRUE) + 0.5)),
# 	TIME = slider(min=0, max=max(df$time))
# )


# 
# ggplot(subset(df, infections <= max(dfv$v, na.rm=TRUE)),
# 			 aes(x=time, y=population, fill=as.factor(infections))) +
# 	scale_fill_grey() +
# 	geom_area(col=NA)

df2 = df %>%
	group_by(time) %>%
	summarise(N = sum(population), I = sum((0:max_i)*population), mean=I/N) %>%
	gather(type, pop, -time) %>%
	arrange(time, type)

ggplot(df2, aes(x=time, y=pop, col=type)) + geom_line()

ggplot(subset(df2, type=="mean"), aes(x=time, y=pop)) + geom_line()



KLD <- function(par, dat) {
  negbin <- dnbinom(x=0:(length(dat)-1), mu=par[1], size=par[2])
  KLD <- sum(dat*pmax(log(dat/negbin), 0))
  return(KLD)
}

library(optimx)
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
  return(optim(par=c(mu=mu.start, size=size.start), fn=KLD, dat=dat, method="L-BFGS-B", lower=c(0,0), upper=(Inf, Inf)))
}

fits = df %>%
	group_by(time) %>% 
	summarize(N = sum(population), I = sum((0:parms$max_i)*population), mean=I/N) %>%
	mutate(mu=NA, size=NA, KLD=NA, convergence=NA)

i = 1
fit = NegBinFit(df[df$time == fits[i,]$time,]$population)
fits[i,]$mu = fit$par[["mu"]]
fits[i,]$size = fit$par[["size"]]
fits[i,]$KLD = fit$value
fits[i,]$convergence = fit$convergence

pbar = txtProgressBar(min=0, max=nrow(fits)-1, initial=1)
for(i in 2:(nrow(fits)-1)) {
	last_fit = max(which(fits$convergence==0))
	fit = try(NegBinFit(df[df$time == fits[i,]$time,]$population,
									mu.start=fits[last_fit,]$mu,
									size.start=fits[last_fit,]$size))
	if (class(fit) == "try-error") {
		print(fit)
		break()
	}
	fits[i,]$mu = fit$par[["mu"]]
	fits[i,]$size = fit$par[["size"]]
	fits[i,]$KLD = fit$value
	fits[i,]$convergence = fit$convergence
	setTxtProgressBar(pbar, i)
}

a = plyr::ddply(df, "time", function(x) {
	v = try(NegBinFit(x$population), silent=TRUE)
	if(class(v) == "try-error") {
		out = data.frame(N = sum(x$population), I = sum(0:100*x$population), mu1=NA, size=NA, value=NA, convergence=NA, var=NA, mu=NA, k=NA) 
	} else {
		out = data.frame(N = sum(x$population), I = sum(0:100*x$population), mu1=v$par[["mu"]], size=v$par[["size"]], value=v$value, convergence=v$convergence, var=var(x$population))
		out$mu = out$I/out$N
		out$k = out$mu^2 / (out$var - out$mu)
	} 
	return(out)
})

aa = a %>% group_by(time) %>% mutate(I_sim = sum(rnbinom(n=N, mu=mu, size=size)))

plot(aa$time, aa$value, type="l")

plot(aa$time, aa$size, type="l", log="y")

plot(aa$time, aa$I, type="l")
lines(aa$time, aa$I_sim, col="red")

b = plyr::ddply(df, "time", function(x) {
	v = try(fitdistr(x$population, "Negative Binomial"))
	if(class(v) == "try-error") {
		out = c(mu=NA, size=NA, value=NA) 
	} else {
		out = c(v$par, value=v$value)
	} 
	return(v)
}, .progress="time")

plot(a$time, a$mu, type="l")
mm = subset(df2, type=="mean")
lines(x=mm$time, y=mm$pop, col="red")

aa = a[complete.cases(a), ]
aa = aa[aa$convergence == 0, ]
plot(a$time, a$size, type="l")
plot(aa$time, aa$size, type="l", col="red")


