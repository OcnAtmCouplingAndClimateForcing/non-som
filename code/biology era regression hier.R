library(tidyverse)
library(zoo)
library(nlme)
library(MARSS)
library(MuMIn)
library(ggpubr)
library(reshape2)
library(glmmTMB)
library(rstanarm)
library(lme4)
library(rstan)
library(ggthemes)
library(tidybayes)
library(bayesplot)

# download PDO / NPGO and process
# download.file("http://jisao.washington.edu/pdo/PDO.latest", "~pdo")
names <- read.table("~pdo", skip=30, nrows=1, as.is = T)
pdo <- read.table("~pdo", skip=32, nrows=119, fill=T, col.names = names)
pdo$YEAR <- 1900:(1899+nrow(pdo)) # drop asterisks!
pdo <- pdo %>%
  gather(month, value, -YEAR) %>%
  arrange(YEAR)

# download.file("http://www.oces.us/npgo/npgo.php", "~npgo")

npgo <- read.table("~npgo", skip=10, nrows=828, fill=T, col.names = c("Year", "month", "value"))

# calculate NDJFM means for each index
pdo$win.yr <- ifelse(pdo$month %in% c("NOV", "DEC"), pdo$YEAR+1, pdo$YEAR)

# limit to winter months only
pdo <- pdo %>%
  filter(month %in% c("NOV", "DEC", "JAN", "FEB", "MAR"))

win.pdo <- tapply(pdo$value, pdo$win.yr, mean)

npgo$win.yr <- ifelse(npgo$month %in% 11:12, npgo$Year+1, npgo$Year)

# limit to winter months only
npgo <- npgo %>%
  filter(month %in% c(11,12,1:3))

win.npgo <- tapply(npgo$value, npgo$win.yr, mean)

# and smoothed (2yr) values of each
win.npgo <- rollapply(win.npgo, 2, mean, align="right", fill=NA)
names(win.npgo) <- 1950:2019
win.pdo <- rollapply(win.pdo, 2, mean, align="right", fill=NA)
names(win.pdo) <- 1900:2019

# load five "non-salmon" data sets: EBS groundfish recruitment, GOA crustaceans/fish,Farallon seabirds, CalCOFI ichthyo,
# CCC seabirds
dat <- read.csv("data/farallon.sbrd.biol.csv", row.names = 1)
# add era term
dat$era <- as.factor(ifelse(dat$year <= 1988, 1, 2))

# and pdo/npgo
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]

# reshape with year, era, and pdo and npgo as the grouping variables
melted <- melt(dat, id.vars = c("year","pdo","era","npgo"))
melted$variable_era = paste0(melted$era,melted$variable)
# standardize all the time series by variable -- so slopes are on same scale
m1 = dplyr::group_by(melted, variable) %>%
  mutate(scale_x = scale(value))
m1$system <- "Central California Current"

####
# CalCOFI
dat <- read.csv("data/calcofi.biol.csv", row.names=1)

# examine distributions
look <- dat %>%
  gather(key, value, -year)

ggplot(look, aes(value)) +
  geom_histogram() +
  facet_wrap(~key, scales="free")

dat$era <- as.factor(ifelse(dat$year <= 1988, 1, 2))

# and pdo/npgo
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]

# reshape with year, era, and pdo and npgo as the grouping variables
melted <- melt(dat, id.vars = c("year","pdo","era","npgo"))
melted$variable_era = paste0(melted$era,melted$variable)
# standardize all the time series by variable -- so slopes are on same scale
m2 = dplyr::group_by(melted, variable) %>%
  mutate(scale_x = scale(value))
m2$system <- "Southern California Current"

#########

dat <- read.csv("data/goa.biol.csv")
colnames(dat)[1] <- "year"

# examine distributions
look <- dat %>%
  gather(key, value, -year)

ggplot(look, aes(value)) +
  geom_histogram() +
  facet_wrap(~key, scales="free")

dat$era <- as.factor(ifelse(dat$year <= 1988, 1, 2))

# and pdo/npgo
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]

# reshape with year, era, and pdo and npgo as the grouping variables
melted <- melt(dat, id.vars = c("year","pdo","era","npgo"))
melted$variable_era = paste0(melted$era,melted$variable)
# standardize all the time series by variable -- so slopes are on same scale
m3 = dplyr::group_by(melted, variable) %>%
  mutate(scale_x = scale(value))
m3$system <- "Gulf of Alaska"

#######
dat <- read.csv("data/ebs.biol.data.csv")
dat[,2:5] <- sqrt(dat[,2:5])

# examine distributions
look <- dat %>%
  gather(key, value, -year)

ggplot(look, aes(value)) +
  geom_histogram() +
  facet_wrap(~key, scales="free")

dat$era <- as.factor(ifelse(dat$year <= 1988, 1, 2))

# and pdo/npgo
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]

# reshape with year, era, and pdo and npgo as the grouping variables
melted <- melt(dat, id.vars = c("year","pdo","era","npgo"))
melted$variable_era = paste0(melted$era,melted$variable)
melted$value <- as.numeric(melted$value)

# standardize all the time series by variable -- so slopes are on same scale
m4 = dplyr::group_by(melted, variable) %>%
  mutate(scale_x = scale(value))
m4$system <- "Bering Sea"

#######
dat <- read.csv("data/ncc.biol.dat.csv")
dat[,2:5] <- sqrt(dat[,2:5])

# examine distributions (only one variable at this point!)
ggplot(dat, aes(value)) +
  geom_histogram() 

dat$era <- as.factor(ifelse(dat$year <= 1988, 1, 2))

# and pdo/npgo
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]

melted <- dat
melted$variable_era = paste0(melted$era,melted$variable)
melted$value <- as.numeric(melted$value)

# standardize all the time series by variable -- so slopes are on same scale
m5 = dplyr::group_by(melted, variable) %>%
  mutate(scale_x = scale(value))
m5$system <- "Northern California Current"

#######
# combine

melted <- rbind(m1, m2, m3, m4, m5)
melted$year <- as.numeric(melted$year)
melted$variable <- as.factor(melted$variable)
melted$variable_era <- as.factor(melted$variable_era)

model.data = data.frame()

levels.syst <- as.factor(unique(melted$system))

for(s in levels.syst) {

  # s <- levels.syst[1]

  temp <- melted %>%
    filter(system==s)
  temp <- na.omit(temp)

  # create data for stan
  temp$variable = as.character(temp$variable)
  temp$variable = as.factor(temp$variable)

  stan_data = list(era = as.numeric(temp$era),
                   y = temp$pdo,
                   variable = as.numeric(temp$variable),
                   n = nrow(temp),
                   n_levels = max(as.numeric(temp$variable)),
                   x = temp$scale_x)

  mod = stan(file="models/mod.stan", data=stan_data, chains=3, warmup=4000, iter=6000,thin=2,
    pars = c("beta","mu_beta","ratio","mu_ratio","sigma_beta","sigma_ratio",
      "exp_mu_ratio","exp_ratio","pred"),
    control=list(adapt_delta=0.99, max_treedepth=20))

  pars = rstan::extract(mod,permuted=TRUE)

  model.data = rbind(model.data,
    data.frame(system=s,
      ratio=100*(pars$mu_ratio)))

  temp$pred = apply(pars$pred,2,mean)

  png(file=paste0("plots/diag/obspred_SI_PDO_",s,".png"))
  g = ggplot(temp, aes(year,pdo,col=era)) + geom_point() +
    facet_wrap(~ variable, scale="free_y") + geom_line(aes(year,pred)) +
    ggtitle(paste0("PDO: Region ",s))
  print(g)
  dev.off()

  png(file=paste0("plots/diag/slopes_SI_PDO_",s,".png"))
  g = ggplot(temp, aes(scale_x,pred,col=era)) + geom_point() +
    facet_wrap(~ variable, scale="free") + xlab("Covariate") +
    ylab("Predicted") + ggtitle(paste0("PDO: Region ",s))
  print(g)
  dev.off()

  # create an array for caterpillar plots by variable
  draws = rstan::extract(mod,permuted=FALSE)
  par_names = dimnames(draws)$parameters
  par_names[grep("^mu_ratio", par_names)] = "global mean"
  par_names[grep("^ratio", par_names)] = levels(temp$variable)
  dimnames(draws)$parameters = par_names
  idx = which(par_names %in% c("global mean",levels(temp$variable)))

  png(paste0(file="plots/SI_PDO_",s,".png"))
  g = mcmc_intervals(draws,
    pars = c("global mean",levels(temp$variable))) +
    geom_vline(xintercept=0,col="red",linetype="dashed") +
    ggtitle(paste0("PDO: Region ",s)) +
    xlab("Avg ratio: Era 2 slope / Era 1 slope") +
    theme_linedraw()
  print(g)
  dev.off()

  # Fit null model to look at autocorrelation between model with and without 2 slopes
  mod0 = stan(file="models/mod0.stan", data=stan_data, chains=3, warmup=4000, iter=6000,thin=2,
             pars = c("beta","mu_beta","sigma_beta","pred"),
             control=list(adapt_delta=0.99, max_treedepth=20))
  pars0 = rstan::extract(mod0,permuted=TRUE)
  resid_df = data.frame("variable"=temp$variable, "year"=temp$year,
                        "x"=temp$scale_x)
  resid_df$pred = apply(pars$pred,2,mean) # predictions from 2slope model
  resid_df$pred0 = apply(pars0$pred,2,mean)# predictions from 1slope model
  resid_df$resid = temp$value - resid_df$pred
  resid_df$resid0 = temp$value - resid_df$pred0
  resids = group_by(resid_df, variable) %>%
    arrange(year) %>%
    summarize(acf2 = unlist(acf(resid,plot=FALSE))[2],
              acf1 = unlist(acf(resid0,plot=FALSE))[2])
  resids$response = "PDO"
  write.csv(resids,file=paste0("output/PDO_",s,"_resid",".csv"))
}


# order the systems north-south
model.data$order <- ifelse(model.data$system=="Bering Sea", 1,
                           ifelse(model.data$system=="Gulf of Alaska", 2,
                                  ifelse(model.data$system=="Northern California Current", 3,
                                         ifelse(model.data$system=="Central California Current", 4, 5))))

model.data$system <- reorder(model.data$system, model.data$order)

pdo.biol.data <- model.data

# save for future reference
write.csv(pdo.biol.data, "output/pdo_biology_model_data.csv")

pdo.biol.data <- read.csv("output/pdo_biology_model_data.csv", row.names=1)

#################
## and the same thing for npgo
model.data <- data.frame()
for(s in levels.syst) {

  # s <- levels.syst[1]

  temp <- melted %>%
    filter(system==s)
  temp <- na.omit(temp)
  # create data for stan
  temp$variable = as.character(temp$variable)
  temp$variable = as.factor(temp$variable)

  stan_data = list(era = as.numeric(temp$era),
                   y = temp$npgo,
                   variable = as.numeric(temp$variable),
                   n = nrow(temp),
                   n_levels = max(as.numeric(temp$variable)),
                   x = temp$scale_x)

  mod = stan(file="models/mod.stan", data=stan_data, chains=3, warmup=4000, iter=6000,thin=2,
    pars = c("beta","mu_beta","ratio","mu_ratio","sigma_beta","sigma_ratio",
      "exp_mu_ratio","exp_ratio","pred"),
    control=list(adapt_delta=0.99, max_treedepth=20))

  pars = rstan::extract(mod,permuted=TRUE)

  model.data = rbind(model.data,
    data.frame(system=s, ratio=100*(pars$mu_ratio)))

  temp$pred = apply(pars$pred,2,mean)

  png(file=paste0("plots/diag/obspred_SI_NPGO_",s,".png"))
  g = ggplot(as.data.frame(temp), aes(year,npgo,col=era)) + geom_point() +
    facet_wrap(~ variable, scale="free_y") + geom_line(aes(year,pred)) +
    ggtitle(paste0("NPGO: Region ",s))
  print(g)
  dev.off()

  png(file=paste0("plots/diag/slopes_SI_NPGO_",s,".png"))
  g = ggplot(temp, aes(scale_x,pred,col=era)) + geom_point() +
    facet_wrap(~ variable, scale="free") + xlab("Covariate") +
    ylab("Predicted") + ggtitle(paste0("NPGO: Region ",s))
  print(g)
  dev.off()

  # create an array for caterpillar plots by variable
  draws = rstan::extract(mod,permuted=FALSE)
  par_names = dimnames(draws)$parameters
  par_names[grep("^mu_ratio", par_names)] = "global mean"
  par_names[grep("^ratio", par_names)] = levels(temp$variable)
  dimnames(draws)$parameters = par_names
  idx = which(par_names %in% c("global mean",levels(temp$variable)))

  png(paste0(file="plots/SI_NPGO_",s,".png"))
  g = mcmc_intervals(draws,
    pars = c("global mean",levels(temp$variable))) +
    geom_vline(xintercept=0,col="red",linetype="dashed") +
    ggtitle(paste0("NPGO: Region ",s)) +
    xlab("Avg ratio: Era 2 slope / Era 1 slope") +
    theme_linedraw()
  print(g)
  dev.off()

  # Fit null model to look at autocorrelation between model with and without 2 slopes
  mod0 = stan(file="models/mod0.stan", data=stan_data, chains=3, warmup=4000, iter=6000,thin=2,
              pars = c("beta","mu_beta","sigma_beta","pred"),
              control=list(adapt_delta=0.99, max_treedepth=20))
  pars0 = rstan::extract(mod0,permuted=TRUE)
  resid_df = data.frame("variable"=temp$variable, "year"=temp$year,
                        "x"=temp$scale_x)
  resid_df$pred = apply(pars$pred,2,mean) # predictions from 2slope model
  resid_df$pred0 = apply(pars0$pred,2,mean)# predictions from 1slope model
  resid_df$resid = temp$value - resid_df$pred
  resid_df$resid0 = temp$value - resid_df$pred0
  resids = group_by(resid_df, variable) %>%
    arrange(year) %>%
    summarize(acf2 = unlist(acf(resid,plot=FALSE))[2],
              acf1 = unlist(acf(resid0,plot=FALSE))[2])
  resids$response = "NPGO"
  write.csv(resids,file=paste0("output/NPGO_",s,"_resid",".csv"))
}


# order the systems north-south
model.data$order <- ifelse(model.data$system=="Bering Sea", 1,
                           ifelse(model.data$system=="Gulf of Alaska", 2,
                                  ifelse(model.data$system=="Northern California Current", 3,
                                         ifelse(model.data$system=="Central California Current", 4, 5))))

model.data$system <- reorder(model.data$system, model.data$order)

npgo.biol.data <- model.data

# save for future reference
write.csv(npgo.biol.data, "output/npgo_biology_model_data.csv")

npgo.biol.data <- read.csv("output/npgo_biology_model_data.csv", row.names=1)
# Caterpillar Plot ===============================
# Helper Functions
q.50 <- function(x) { return(quantile(x, probs=c(0.25,0.75))) }
q.90 <- function(x) { return(quantile(x, probs=c(0.05,0.95))) }
q.95 <- function(x) { return(quantile(x, probs=c(0.025,0.975))) }

head(npgo.biol.data)
head(pdo.biol.data)

# Combine dataframes
npgo.biol.data$var <- "NPGO"
pdo.biol.data$var <- "PDO"

# now add in salmon results
salmon <- read.csv("output/Salmon/list.mu.ratio.csv", row.names=1)

unique(salmon$region)

salmon$region <- as.factor(ifelse(salmon$region=="EBS", "Bering Sea",
                        ifelse(salmon$region=="GOA", "Gulf of Alaska", "Northern California Current")))

salmon$var <- as.factor(ifelse(salmon$var=="pdo2a", "PDO", "NPGO"))

salmon$order <- as.factor(ifelse(salmon$region=="Bering Sea", 1, 
                                 ifelse(salmon$region=="Gulf of Alaska", 2, 3)))

salmon <- salmon %>%
  select(region, value, order, var)

names(salmon)[c(1,2)] <- c("system", "ratio")
salmon$var.order <- ifelse(salmon$var=="PDO", 1, 2)


# all.data <- rbind(pdo.biol.data, npgo.biol.data, salmon)
all.data <- rbind(pdo.biol.data, npgo.biol.data)
all.data$var.order <- ifelse(all.data$var=="PDO", 1, 2)
all.data$var <- reorder(all.data$var, all.data$var.order)

# order the systems north-south
all.data$order <- ifelse(all.data$system=="Bering Sea", 1,
                           ifelse(all.data$system=="Gulf of Alaska", 2,
                                  ifelse(all.data$system=="Northern California Current", 3,
                                         ifelse(all.data$system=="Central California Current", 4, 5))))

all.data$system <- reorder(all.data$system, all.data$order)
# all.data <- rbind(all.data, salmon) # uncomment to combine DFs


# colorblind...
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

biol.plt <- ggplot(all.data, aes(x=reorder(system, desc(system)), y=ratio, fill=system)) +
  theme_linedraw() +
  scale_fill_manual(values=cb[c(6,3,4,2,8)],
                    labels=c("Bering Sea", "Gulf of Alaska",
                             "Northern Cal. Curr.", "Central Cal. Curr.", "Southern Cal. Curr.")) +
  # scale_fill_colorblind() +
  # scale_fill_tableau() +
  # scale_fill_brewer(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
  # geom_eye() +

  geom_violin(alpha = 0.75, lwd=0.1, scale='width') +
  stat_summary(fun.y="q.90", colour="black", geom="line", lwd=1) +
  stat_summary(fun.y="median", colour="black", size=2, geom="point", pch=21) +
  facet_wrap(~var, ncol=1) +
  ylab("Ratio: Era 1 slope / Era 2 slope (%)") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_line(size=0),
        legend.title = element_blank(), legend.position = c(0.15,0.15)) +
  geom_hline(aes(yintercept=1), color="red", linetype="dotted", size=1) +
  coord_flip(ylim=c())

biol.plt

salmon.plt <- ggplot(salmon, aes(x=reorder(system, desc(system)), y=ratio, fill=system)) +
  theme_linedraw() +
  scale_fill_manual(values=cb[c(6,3,4,2,8)],
                    labels=c("Bering Sea", "Gulf of Alaska",
                             "Northern Cal. Curr.", "Central Cal. Curr.", "Southern Cal. Curr.")) +
  # scale_fill_colorblind() +
  # scale_fill_tableau() +
  # scale_fill_brewer(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
  # geom_eye() +
  
  geom_violin(alpha = 0.75, lwd=0.1, scale='width') +
  stat_summary(fun.y="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun.y="q.50", colour="black", geom="line", lwd=1.5) +
  stat_summary(fun.y="median", colour="black", size=2, geom="point", pch=21) +
  facet_wrap(~var, ncol=1) +
  ylab("Log ratio: Era 1 slope / Era 2 slope") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_line(size=0),
        legend.title = element_blank(), legend.position = c(0.15,0.15)) +
  geom_hline(aes(yintercept=1), color="red", linetype="dotted", size=1) +
  coord_flip(ylim=c())

ggsave("plots/salmon plot.png")

cat.plt <- ggplot(all.data, aes(x=reorder(system, desc(order)), y=ratio/100, fill=system)) +
             theme_linedraw() +
             # scale_fill_colorblind() +
             # scale_fill_tableau() +
             # scale_fill_brewer(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
             # geom_eye() +
             scale_fill_manual(values=cb[c(6,3,4,2,8)],
                    labels=c("Bering Sea", "Gulf of Alaska",
                             "Northern Cal. Curr.", "Central Cal. Curr.", "Southern Cal. Curr.")) +
             geom_violin(alpha = 0.75, lwd=0.1, scale='width') +
             stat_summary(fun.y="q.95", colour="black", geom="line", lwd=0.75) +
             stat_summary(fun.y="q.50", colour="black", geom="line", lwd=1.5) +
             stat_summary(fun.y="median", colour="black", size=2, geom="point", pch=21) +
             facet_wrap(~var, ncol=1) +
             ylab("Avg ratio: Era 1 slope / Era 2 slope") +
             theme(axis.text.y = element_blank(), axis.title.y = element_blank(), legend.title = element_blank()) +
             geom_hline(aes(yintercept=1), color="red", linetype="dotted", size=1) +
             coord_flip(ylim=c(-3,3))

cat.plt

ggsave("biol regression change pdo-npgo slope_cater.png", plot=cat.plt,
         height=7, width=7, units="in", dpi=300)

# Separate by

cat.plt.pdo <- all.data %>% filter(var=='PDO') %>% ggplot(aes(x=system, y=ratio/100, fill=system)) +
  theme_linedraw() +
  scale_fill_tableau() +
  # geom_eye() +
  geom_violin(alpha = 0.75, lwd=0.1, scale='width') +
  stat_summary(fun.y="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun.y="q.50", colour="black", geom="line", lwd=1.5) +
  stat_summary(fun.y="median", colour="black", size=2, geom="point", pch=21) +
  facet_wrap(~var, ncol=1) +
  ylab("Avg ratio: Era 1 slope / Era 2 slope") +
  theme(axis.text.y = element_blank(), legend.position='top') +
  geom_hline(aes(yintercept=1), color="red", linetype="dotted", size=1) +
  coord_flip(ylim=c(-3,3))

cat.plt.npgo <- all.data %>% filter(var=='NPGO') %>% ggplot(aes(x=system, y=ratio/100, fill=system)) +
  theme_linedraw() +
  scale_fill_tableau() +
  # geom_eye() +
  geom_violin(alpha = 0.75, lwd=0.1, scale='width') +
  stat_summary(fun.y="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun.y="q.50", colour="black", geom="line", lwd=1.5) +
  stat_summary(fun.y="median", colour="black", size=2, geom="point", pch=21) +
  facet_wrap(~var, ncol=1) +
  ylab("Avg ratio: Era 1 slope / Era 2 slope") +
  theme(axis.text.y = element_blank(), legend.position="none") +
  geom_hline(aes(yintercept=1), color="red", linetype="dotted", size=1) +
  coord_flip(ylim=c(-3,3))

# Plot Combined with Sepearte
cat.plt.2 <- plot_grid(cat.plt.pdo, cat.plt.npgo, ncol=1, rel_heights = c(1.1,1))
ggsave("biol regression change pdo-npgo slope_cater2.png", plot=cat.plt.2,
       height=7, width=7, units="in", dpi=300)

#Plot: Facet by System =====================
cat.plt.3 <- ggplot(all.data, aes(x=var, y=ratio/100, fill=var)) +
  theme_linedraw() +
  # scale_fill_colorblind() +
  scale_fill_tableau() +
  # geom_eye() +

  geom_violin(alpha = 0.75, lwd=0.1, scale='width') +
  stat_summary(fun.y="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun.y="q.50", colour="black", geom="line", lwd=1.5) +
  stat_summary(fun.y="median", colour="black", size=2, geom="point", pch=21) +
  facet_wrap(~system, ncol=1, scales='free_y') +
  ylab("Avg ratio: Era 1 slope / Era 2 slope") +
  xlab("") +
  theme(axis.text.y = element_blank(),
        legend.position = 'top') +
  geom_hline(aes(yintercept=1), color="red", linetype="dotted", size=1) +
  coord_flip(ylim=c(-3,3))


cat.plt.3
ggsave("biol regression change pdo-npgo slope_cater3.png", plot=cat.plt.3,
       height=7, width=7, units="in", dpi=300)


