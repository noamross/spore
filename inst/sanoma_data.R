library(readr)
library(dplyr)
library(lubridate)
library(magrittr)
library(tidyr)
library(ggplot2)
library(ggmap)
library(stringi)
library(rgdal)
library(cmp)

stems = read.csv('data/host.stems_sonoma.csv')
stems %<>% tbl_df %>%
  mutate(date = dmy(date), year = year(date), plotname = stri_extract_first_regex(plotid, "[A-Z]+"),
            plotnum = stri_extract_first_regex(plotid, "\\d+"))


canker_counts = stems %>%
  filter(status=="Alive") %>%
  group_by(date, tagnumber) %>%
  summarize(ncanker = sum(canker), plotid = plotid[1], plotname = plotname[1], plotnum = plotnum[1], species=species[1], symplfct=symplfct[1], lidelfct = lidelfct[1], year=year[1]) %>%
  group_by()

canker_total = canker_counts %>%
  filter(year >= 2010) %>%
  group_by(year, species) %>%
  summarize(infected = sum(ncanker > 0), meancank = mean(ncanker), varcank = var(ncanker), meanvar = varcank/meancank, trees = n(), fracinf=infected/trees) %>%
  group_by() %>%
  arrange(species, year)

canker_stats = canker_counts %>%
  filter(year >= 2010) %>%
  group_by(plotid, year, species) %>%
  summarize(infected = sum(ncanker > 0), meancank = mean(ncanker), varcank = var(ncanker), meanvar = varcank/meancank, trees = n(), fracinf=infected/trees) %>%
  group_by() %>%
  arrange(meanvar)



ggplot(filter(canker_stats, species != "UMCA", infected > 0), aes(x=trees, y =log(meanvar)), col=species) + geom_point()
ggplot(filter(canker_stats, infected > 0), aes(x=meanvar, fill=as.factor(year))) + geom_histogram(binwidth=0.25, col="black") + facet_wrap(~year, ncol=1) + scale_x_continuous(breaks=0:10)


leaf_stats = canker_counts %>%
  filter(species == "UMCA", symplfct != -9999) %>%
  group_by(plotid, year) %>%
  summarize(trees = n(), infected = sum(symplfct > 0), meanlf = mean(symplfct),
            varlf = var(symplfct), mvar = meanlf/varlf) %>%
  group_by() %>%
  arrange(desc(infected))

ggplot(filter(leaf_stats, year==2012), aes(x=meanlf, y =log(mvar))) + geom_point()


lide_stats = canker_counts %>%
  filter(species == "LIDE", lidelfct != -9999) %>%
  group_by(plotid, year) %>%
  summarize(infected = sum(lidelfct > 0), meanlf = mean(lidelfct),
            varlf = var(lidelfct), mvar = meanlf/varlf)

# ggplot(filter(lide_stats, year %in% 2010:2012), aes(x=meanlf, y =log(mvar))) + geom_point()
#
# ggplot(filter(canker_counts, year==2012, species == "LIDE", plotid != "HOOD03"), aes(x = lidelfct, fill=plotid)) + geom_histogram() + facet_wrap(~plotid)
#
# ggplot(filter(canker_counts, year==2012, species == "UMCA"), aes(x = symplfct, fill=plotid)) + geom_histogram() + facet_wrap(~plotid)
#

locations = stems %>%
  filter(!is.na(x), !is.na(y)) %>%
  group_by(plotid) %>%
  summarize(x = x[1], y=y[1]) %>%
  mutate(plotname = stri_extract_first_regex(plotid, "[A-Z]+"),
            plotnum = stri_extract_first_regex(plotid, "\\d+"))

sputm <- SpatialPoints(data.frame(locations$x, locations$y), proj4string=CRS("+proj=utm +zone=10 +datum=WGS84"))
spgeo <- as.data.frame(spTransform(sputm, CRS("+proj=longlat +datum=WGS84")))
names(spgeo) <- c("lon", "lat")
locations = cbind(locations, spgeo)
bb = make_bbox(locations$lon, locations$lat)
map = get_map(bb)
ggmap(map) + geom_point(data = locations, mapping = aes(x = lon, y = lat, fill = plotname), col = "black", size = 3, pch=21) + theme(legend.position = "none")


# These are 15 x 15 m plots
# See http://dx.doi.org/10.1111/j.1365-2745.2008.01376.x

library(cmp)
library(rethinking)
fit_data = canker_counts %>%
  filter(year == 2011, species == "LIDE", plotname=="JOHNS") %>%
  as.data.frame

dnbinom2 = function(x, mu, size, log = FALSE) {dnbinom(x=x, mu=mu, size=size, log=log)}

sano = canker_counts %>%
  filter(year >= 2010, species != "UMCA") %>%
  group_by(year, species, plotname) %>%
  summarize(cankered = sum(ncanker > 0), meancank = mean(ncanker), meaninfcank = mean(ncanker[ncanker > 0]), varcank = var(ncanker), meanvar = varcank/meancank, trees = n(), fraccank=cankered/trees) %>%
  group_by() %>%
  arrange(meanvar)

counter = function(x) {
  tab = table(x)
  out = cbind(as.numeric(names(tab)), tab)
  rownames(out) <- NULL
  colnames(out) <- NULL
  return(out)
  }


cmp_fits = canker_counts %>%
  filter(year >= 2010) %>%
  group_by(year, plotname, species) %>%
  do({
    counts = counter(.$ncanker)
    fit_cmp = cmp_fit(counts)
    fit_pois = pois_fit(counts)
    fit_nb = nb_fit(counts)
    data.frame(cankered = sum(.$ncanker > 0), meancank = mean(.$ncanker), meaninfcank = mean(.$ncanker[.$ncanker > 0]), varcank = var(.$ncanker), meanvar = var(.$ncanker)/mean(.$ncanker), trees = length(.$ncanker), fraccank=sum(.$ncanker > 0)/length(.$ncanker),
      lambda_cmp = fit_cmp$lambda, nu_cmp = fit_cmp$nu, z_cmp = fit_cmp$z, ll_cmp = fit_cmp$log.likelihood, lambda_pois = fit_pois$lambda, ll_pois = fit_pois$log.likelihood,
      mu_nb = fit_nb$mu, size_nb = fit_nb$size, ll_nb = fit_nb$log.likelihood)
  }) %>%
  group_by() %>%
  arrange(nu_cmp)


sanoma_model_cmp = alist(
  ncanker ~ dcmp(lambda, nu),
  log(lambda) <- a,
  log(nu) <- b,
  a ~ dnorm(0,10),
  b ~ dnorm(0,10))

sanoma_model_pois = alist(
  ncanker ~ dpois(lambda),
  log(lambda) <- a,
  a ~ dnorm(0,10))

sanoma_model_nb = alist(
  ncanker ~ dnbinom2(mu, size),
  log(mu) <- a,
  log(size) <- b,
  a ~ dnorm(0,10),
  b ~ dnorm(0,10))

fit_data = data.frame(ncanker=sample(c(0,1), 30, replace=TRUE))

fit_sanoma_cmp = map(sanoma_model_cmp, data = fit_data)
fit_sanoma_pois = map(sanoma_model_pois, data = fit_data)
fit_sanoma_nb = map(sanoma_model_nb, data = fit_data)
precis(fit_sanoma_cmp)
precis(fit_sanoma_pois)
precis(fit_sanoma_nb)
compare(fit_sanoma_cmp, fit_sanoma_pois, fit_sanoma_nb, WAIC=TRUE)
AIC(fit_sanoma_cmp, fit_sanoma_pois, fit_sanoma_nb)


counter = function(x) {
  tab = table(x)
  out = cbind(as.numeric(names(tab)), tab)
  rownames(out) <- NULL
  colnames(out) <- NULL
  return(out)
  }
