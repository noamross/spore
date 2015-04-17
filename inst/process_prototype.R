library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(rlist)
library(noamtools)
library(manipulate)
library(compoisson)
library(scales)

source("R/negbinfit.R")

fname="out_ad.txt"
#con = file(fname)
#out1 = lapply(strsplit(readLines(con), " "), as.numeric)
n_col <- max(count.fields(fname, sep = " "))
output <- read.table(fname ,sep=" ",fill=TRUE, col.names = 1:n_col)
colnames(output) = c("time", as.character(0:(n_col-2)))
output[is.na(output)] = 0
output = output[, colSums(output) !=0]
max_i = ncol(output) - 2

df = output %>%
	gather(infections, population, -time, convert=TRUE) %>%
	arrange(time, infections)

fits = df %>%
	group_by(time) %>% 
	summarize(N = sum(population),
						P = sum((0:max_i)*population),
						S = sum(population[infections == 0]),
						I = N - S,
						mean=P/N,
						variance=var(population),
						variance_mean_ratio = variance/mean)

samp_times = fits$time[findInterval(seq(min(fits$time), max(fits$time), length.out=500), fits$time)]
fits %<>% filter(time %in% samp_times)
df %<>% filter(time %in% samp_times)

rm(output)

fits2 = fits %>% gather(variable, value, -time)
ggplot(subset(fits2, variable %in% c("S", "I", "N", "P")), aes(x=time, y=value, col=variable)) + geom_line() + theme_nr 
ggplot(subset(fits2, variable %in% c("mean")), aes(x=time, y=value, col=variable)) + geom_line() + theme_nr 



storage = new.env()
storage$pars = c(1,1)

#pbar = txtProgressBar(min=1, max=nrow(fits), initial=1)
dfits = df %>%
#	filter(time %in% samp_times) %>%
	group_by(time) %>%
	do({
		fit = com.fit(cbind(.$infections, .$population), .GlobalEnv$storage$pars, maxit=10)
		if(!any(c(fit$lambda, fit$nu) > 1e10 | c(fit$lambda, fit$nu) < 1e-10)) {
			.GlobalEnv$storage$pars = c(fit$lambda, fit$nu)
		}
		fit2 = pois.fit(cbind(.$infections, .$population))
		fit3 = nb.fit(cbind(.$infections, .$population))
		return(data.frame(list(list.filter(fit, .name != "fitted.values"),
													 pois_lambda=fit2$lambda, pois_loglik=fit2$log.likelihood,
													 nb_mu = fit3$mu, nb_size = fit3$size, nb_loglik=fit3$log.likelihood)))
			})

dfits = dfits %>%	mutate(fit_mean = com_mean(lambda, nu, z=z),
				 fit_var = com_var(lambda, nu, z=z),
				 com_AIC = 4 - 2*(log.likelihood),
				 pois_AIC = 2 - 2*pois_loglik,
				 NB_AIC = 4 - 2*nb_loglik,
				 NB_com_AIC_diff = NB_AIC - com_AIC,
				 pois_com_AIC_diff = pois_AIC - com_AIC)

fits = inner_join(fits, dfits)

fitcounts = fits %>%
	group_by(time) %>%
	do(data_frame(infections = 0:max_i,
								CMP_fit = .$N * dcom(0:max_i, .$lambda, .$nu),
								Poisson_fit = .$N * dpois(0:max_i, .$pois_lambda),
								NB_fit = .$N * dnbinom(0:max_i, mu=.$nb_mu, size=.$nb_size))) %>%
	inner_join(df) %>%
 	gather(origin, population, population, CMP_fit, Poisson_fit, NB_fit) 

fits2 = fits %>% gather(variable, value, -time)

ggplot(subset(fits2, variable %in% c("S", "I", "N", "P")), aes(x=time, y=value, col=variable)) + geom_line() + theme_nr + ylim(0,100) 

ggplot(subset(fits2, variable %in% c("nu") & time > 1 & time < 90), aes(x=time, y=value, col=variable)) + geom_line() + theme_nr + geom_hline(yintercept = 1) 

ggplot(subset(fits2, variable %in% c("nb_size")), aes(x=time, y=value, col=variable)) + geom_line() + theme_nr + geom_hline(yintercept = 1)

ggplot(subset(fits2, variable %in% c("variance_mean_ratio") & time > 25), aes(x=time, y=value, col=variable)) +
	geom_line() + 
#	scale_y_log10() +
	theme_nr #+

ggplot(fits[-(1:10),], aes(x=variance_mean_ratio, y=nu)) + geom_point() + theme_nr + scale_x_log10() + scale_y_log10()

ggplot(subset(fits2, variable %in% c("log.likelihood", "pois_loglik", "nb_loglik")), aes(x=time, y=value, col=variable)) + geom_line() + theme_nr

ggplot(subset(fits2, variable %in% c("NB_com_AIC_diff", "pois_com_AIC_diff")), aes(x=time, y=value, col=variable)) + geom_line() + theme_nr 

ggplot(subset(fits2, variable %in% c("com_AIC", "pois_AIC", "NB_AIC") & !is.na(value)), aes(x=time, y=value, col=variable)) + 
	geom_line() + 
	theme_nr #+

ggplot(fits, aes(x=mean, y=lambda)) + geom_path() + geom_abline(intercept=0, slope=1)
ggplot(fits, aes(x=variance, y=fit_var)) + geom_line() + scale_x_log10() + scale_y_log10()

ggplot(subset(fits2, variable %in% c("variance", "fit_var") & !is.na(value)), aes(x=time, y=value, col=variable)) + 
	geom_line() + 
	theme_nr #+


ggplot(subset(fits2, variable %in% c("P", "P_sim", "S", "S_sim", "I", "I_sim")), aes(x=time, y=value, col=variable)) + 
	geom_line() + 
	theme_nr #+
	#scale_y_log10()

manipulate({
	ddf = filter(fitcounts, time == max(fitcounts$time[fitcounts$time <= TIME]));
	ggplot() +
		geom_bar(mapping = aes(x = infections, y=population),
						 data = filter(ddf, origin=="population"),
						 position="identity", stat="identity", col=NA, fill="grey") +
		geom_step(mapping = aes(x = infections - 0.5, y=population, col=origin, linetype=origin),
							data = filter(ddf, origin!="population"), direction="hv") +
		#scale_y_continuous(limits = c(0, 100), oob=rescale_none) +
		scale_x_discrete(limits = c(0, max(ddf$infections[which(ddf$infections > 0)])))
	},
	TIME = slider(min = 0, max = max(df$time))
)

manipulate({
	ddf = filter(fitcounts, time == max(fitcounts$time[fitcounts$time <= TIME]));
	ggplot(ddf, aes(x=infections, y=population, fill=origin)) +
		geom_bar(stat="identity", position="identity", alpha=0.5) +
		scale_x_discrete(limits = c(0, max(ddf$infections[which(ddf$infections > 0)])))
	},
	TIME = slider(min = 0, max = max(df$time))
)

ddf = subset(df, time == max(df$time[df$time <= TIME]))
com.fit(cbind(ddf$infections, ddf$population))

library(noamtools)
proftable("out.prof")
				 