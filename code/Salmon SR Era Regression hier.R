#==================================================================================================
#Project Name: FATE Working Group - Hierarchcial Salmon Breakpoint Regression
#Date: 6.12.19
#
#Purpose: Fit hierarchical stock-recruitment model to salmon data that estimates a breakpoint in PDO and NPGO effects
#  a) Fit to each species separately
#  b) Estimate PDO/NPGO effects pre/post breakpoint
#  c) Population level effects assumed to arise from common distribution across locations, but with
#       region-specific mean and standard deviation.
#
#
#==================================================================================================
#NOTES:
#
#==================================================================================================
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
library(brms)

# Define Workflow Paths ============================================
# *Assumes you are working from the Sergent_Streamflow R project
wd <- getwd()
dir.output <- file.path(wd,"output")
dir.figs <- file.path(wd,"plots")
dir.data <- file.path(wd,"data")
dir.mods <- file.path(wd, "models")


# CONTROL ==========================================================
fit <- TRUE

# MCMC Parameters
n.chains <- 3
n.iter <- 5e3
n.thin <- 5

# Select Species 
fit.species <- c("Sockeye","Pink","Chum")[1]

# Whether to fit a model with PDO or NPGO 
vars <- c("pdo1","pdo2a","pdo2b",
          "pdo3","npgo","npgo2a",
          "npgo2b","npgo3")

var <- vars[1]



# Define Output Folders ============================================
dir.output <- file.path(dir.output, paste0(fit.species,"-",var))
dir.create(dir.output)



# Load/Compile Data ================================================
# Stock-recruitment data
dat <- read.csv(file.path(dir.data, "salmon run dat.csv"), header=TRUE, stringsAsFactors=FALSE)
dat.sr <- dat %>% filter(species==fit.species)
dat.sr$ln.rps <- log(dat.sr$recruits/dat.sr$spawners)

# Environmental data
dat.env <- read.csv(file.path(dir.data, "winter pdo npgo various smoothings.csv"), header=TRUE, stringsAsFactors=FALSE)

dat.env.fit <- dat.env %>% select(year, var)

# Bind Environmental data
dat.input <- dat.sr %>% left_join(dat.env.fit, by=c("entry.yr"="year"))
# Rename environmental variable
names(dat.input)[ncol(dat.input)] <- "env.var"

# Remove any missing environmental observations
dat.input <- dat.input %>% filter(!is.na(ln.rps) & !is.na(spawners) & !is.nan(ln.rps) & !is.nan(spawners))

# Metadata
regions <- unique(dat.input$region)
n.regions <- length(regions)


# Fit brms model ===================================================
# fit.brm <- brm(ln.rps ~ (1|stock) + spawners:stock + env.var:era + (env.var:era|region/stock), 
#                data=dat.input,
#                iter=1000, thin=2,
#                chains=3,
#                family=gaussian(),
#                # prior=prior(normal(0,1), class=b))
#                # prior=c(prior(normal(0,1), class=b),
#                #         # prior(normal(0,10), class=Intercept),
#                #         prior(normal(0,10), class=sd),
#                #         # prior(normal(0,100), class=sds),
#                #         prior(normal(0,10), class=sds),
#                #         prior(normal(0,10), class=sigma)),
#                sample_prior=TRUE,
#                control = list(adapt_delta = 0.99))
# 
# summary(fit.brm)

# Create Input Data ================================================
stocks <- unique(dat.input$stock)
n.stocks <- length(stocks)

#Required Attributs of the species
stock.regions <- unique(dat.input$region)
n.stock.regions <- length(stock.regions)

stock.years <- min(dat.input$brood.yr):max(dat.input$brood.yr)
n.stock.years <- length(stock.years)


#Create Data Objects
maxN <- n.stock.years
S <- n.stocks
N <- vector(length=n.stocks)
R <- n.stock.regions #Number of Regions
region <- vector(length=n.stocks)
# K <- 2 #Number of covariates PDO, NPGO
# covars <- array(dim=c(n.stocks,100,K)) #We will start with a temporary length of 100 years then trim down to the max N
covar <- array(data=0, dim=c(n.stocks,maxN))
era <- array(data=0, dim=c(n.stocks,maxN))

#Year Pointer - For Referencing Coefficient Values
years <- array(data=0, dim=c(n.stocks,maxN))

#Ricker Parameters
ln_rps <- array(data=0, dim=c(n.stocks, maxN))
spawn <- array(data=0, dim=c(n.stocks, maxN))

p <- 1
for(p in 1:n.stocks) {
  
  #Retreive Data =================
  temp.stock <- stocks[p]
  dat.temp <- dat.input %>% filter(stock==temp.stock) %>% arrange(brood.yr)
  
  #Assign STAN Inputs ============
  N[p] <- nrow(dat.temp) # Number of years in stock-recruitment time series for stock
  region[p] <- which(stock.regions==unique(dat.temp$region)) # Pointer to region for stock
  
  ln_rps[p, 1:N[p]] <- dat.temp$ln.rps # Response variable for stock
  spawn[p, 1:N[p]] <- dat.temp$spawn # Spawning abundance for density-dependent effect on stock
  
  years[p,1:N[p]] <- which(stock.years %in% dat.temp$brood.yr ) # Stock-specific brood years
  
  #Assign Covars and Eras ===============
  covar[p,1:N[p]] <- dat.temp$env.var
  era[p,1:N[p]] <- dat.temp$era

}#next p

#Determine maximum length of covariates =====================
# maxN <- max(N)
# temp.regions <- regions[unique(region)]

# Fit Stan Model ===================================================
#Fit the model
if(fit==TRUE) {
  stan.fit <- stan(file=file.path(dir.mods,"hier-Ricker.stan"),
             model_name="Hierarchical-Breakpoint-Ricker",
             data=list("ln_rps"=ln_rps, "spawn"=spawn,
                       "N"=N, "maxN"=maxN,
                       "S"=S, "R"=R,
                       "region"=region,
                       "covar"=covar, "era"=era
                       ),
             chains=n.chains, iter=n.iter, thin=n.thin,
             cores=n.chains, verbose=FALSE,
             seed=101)
  #Save Output
  saveRDS(fit, file=file.path(dir.output,"stan.fit.rds"))
}else {
  stan.fit <- readRDS(file=file.path(dir.output,"stan.fit.rds"))
}

# Plot Output ======================================================
# Someone can continue here.

pars <- extract(stan.fit)

mu_ratios <- data.frame(pars$mu_ratio)
names(mu_ratios) <- regions
exp_mu_ratios <- exp(mu_ratios)

list.exp_mu_ratios <- melt(exp_mu_ratios)

g <- list.exp_mu_ratios %>% ggplot(aes(value, fill=variable)) +
       scale_fill_colorblind() +
       geom_density(alpha=0.5)
g

g2 <- list.exp_mu_ratios %>% ggplot(aes(x=variable, y=value, fill=variable)) +
        scale_fill_colorblind() +
        geom_eye(alpha=0.5) +
        coord_flip() +
        theme(legend.position = 'none')
g2




































