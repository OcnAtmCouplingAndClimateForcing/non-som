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
library(cowplot)
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

melted <- m3
melted$year <- as.numeric(melted$year)
melted$variable <- as.factor(melted$variable)
melted$variable_era <- as.factor(melted$variable_era)

levels.syst <- as.factor(unique(melted$system))
levels.vars <- as.factor(unique(melted$variable)) 

model.preds <- data.frame()
coefs <- data.frame()

for(s in levels.syst) {
  
  # s <- levels.syst[1]
  
  # fitting a separate model to each biology time series
  for(v in levels.vars){ 
    
  temp <- melted %>%
    filter(system==s, variable==v)
  temp <- na.omit(temp)
  
  mod = lm(pdo ~ -1 + scale_x:era, data=temp)
  coefs = rbind(coefs, data.frame(b1 = coef(mod)[1],
                                  b2 = coef(mod)[2]))
  model.preds <- rbind(model.preds,
                       data.frame(year=temp$year, variable=temp$variable, 
                                  obs=temp$pdo,
                                  pred=predict(mod)))
  }}

coefs$variable = levels.vars
  
model.preds$era <- as.factor(ifelse(model.preds$year <=1988, 1, 2))
model.preds$variable <- as.character((model.preds$variable))

png(file=paste0("plots/diag/lm_comparison_obspred_SI_PDO_",s,".png"))
  g = ggplot(model.preds, aes(year,obs,col=era)) + geom_point() +
    facet_wrap(~ variable, scale="free_y") + geom_line(aes(year,pred)) +
    ggtitle(paste0("PDO: Region ",s))
  print(g)
  dev.off()


 