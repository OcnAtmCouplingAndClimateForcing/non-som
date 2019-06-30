#==================================================================================================
#Project Name: FATE Working Group - Plotting and Compiling Output
#Date: 6.25.19
#
#Purpose: Compile output objects with values for 
#
#
#==================================================================================================
#NOTES:
# 
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
library(shinystan)

# Define Workflow Paths ============================================
# *Assumes you are working from the Sergent_Streamflow R project
wd <- getwd()
# dir.output <- file.path(wd,"output")
dir.output <- file.path(wd,"output","CurryTesting")
dir.figs <- file.path(wd,"plots", "Salmon")
dir.data <- file.path(wd,"data")
dir.mods <- file.path(wd, "models")

dir.create(dir.figs)

# CONTROL ==========================================================
read <- TRUE #Whether to read in all model files


# MCMC Parameters
# n.chains <- 3
# n.iter <- 5e4
# n.thin <- 10

# Select Species 
species <- c("Sockeye","Pink","Chum")
n.species <- length(species)


# Whether to fit a model with PDO or NPGO 
vars <- c("pdo2a", "npgo2a")
n.vars <- length(vars)

regions <- c("South","GOA","EBS")
n.regions <- length(regions)

# Create Output Objects ==========================================
list.rhat <- matrix(nrow=0, ncol=3) # Convergence diag
list.neff <- matrix(nrow=0, ncol=3) # Effective Sample Size
list.phi <- matrix(nrow=0, ncol=3) # Autoregressive Coefficient

# Effects Ratio (in log space)
list.mu.ratio <- matrix(nrow=0, ncol=4)
list.sigma.ratio <- matrix(nrow=0, ncol=4)

# Load A Sample Dataset ==========================================
if(read==TRUE) {
s <- 1
for(s in 1:n.species) {
  print(paste("s:",s,"of",n.species))
  temp.species <- species[s]
  
  v <- 1
  for(v in 1:n.vars) {
    print(paste("v:",v,"of",n.vars))
    temp.var <- vars[v]
    
    # Load data
    temp.fit <- readRDS(file=file.path(dir.output, paste0(temp.species,"-",temp.var,".rds")))
    temp.sum <- summary(temp.fit)
    
    temp.pars <- extract(temp.fit)
    # Extract Convergence Diagnostics ==========================
    rhat <- temp.sum$summary[,ncol(temp.sum$summary)]
    neff <- temp.sum$summary[,ncol(temp.sum$summary)-1]
    
    list.rhat <- rbind(list.rhat, data.frame(temp.species, temp.var, rhat))
    list.neff <- rbind(list.neff, data.frame(temp.species, temp.var, neff))
    
    # Extract Autoregression Coeff ===========================================
    list.phi <- rbind(list.phi, data.frame(temp.species, temp.var, temp.pars$phi))
    
    # Extract Regional Ratio Information =====================================
    temp.mu.ratio <- data.frame(temp.pars$mu_ratio)
    names(temp.mu.ratio) <- regions
    temp.mu.ratio.2 <- melt(temp.mu.ratio)
    
    list.mu.ratio <- rbind(list.mu.ratio, data.frame(temp.species, temp.var, 
                                                       temp.mu.ratio.2))
    
    temp.sigma.ratio <- data.frame(temp.pars$sigma_ratio)
    names(temp.sigma.ratio) <- regions
    temp.sigma.ratio.2 <- melt(temp.sigma.ratio)
    
    list.sigma.ratio <- rbind(list.sigma.ratio, data.frame(temp.species, temp.var, 
                                                           temp.sigma.ratio.2))
    
    
  }
}
  
# Save Extracted Objects ==========================================
names(list.rhat) <- c("species","var","value")
names(list.neff) <- c("species","var","value")
names(list.phi) <- c("species","var","value")

write.csv(list.rhat, file=file.path(wd,"output","Salmon","list.rhat.csv")) 
write.csv(list.neff, file=file.path(wd,"output","Salmon","list.neff.csv")) 
write.csv(list.phi, file=file.path(wd,"output","Salmon","list.phi.csv")) 


names(list.mu.ratio) <- c("species","var","region","value")
names(list.sigma.ratio) <- c("species","var","region","value")

write.csv(list.mu.ratio, file=file.path(wd,"output","Salmon","list.mu.ratio.csv")) 
write.csv(list.sigma.ratio, file=file.path(wd,"output","Salmon","list.sigma.ratio.csv")) 

}else {
  list.rhat <- read.csv(list.rhat, file=file.path(wd,"output","Salmon","list.rhat.csv")) 
  list.neff <- read.csv(file=file.path(wd,"output","Salmon","list.neff.csv")) 
  list.phi <- read.csv(file=file.path(wd,"output","Salmon","list.phi.csv")) 
  
  list.mu.ratio <- read.csv(file=file.path(wd,"output","Salmon","list.mu.ratio.csv")) 
  list.sigma.ratio <- read.csv(file=file.path(wd,"output","Salmon","list.sigma.ratio.csv")) 
}
  

# Explore Models with Shiny Stan =======================================
shiny.fit <- readRDS(file=file.path(dir.output, paste0("Chum-npgo2a.rds")))

# launch_shinystan(shiny.fit) #Required internet
# deploy_shinystan(as.shinystan(shiny.fit))
  
# Plot: rhat ===========================================================
g.rhat <- ggplot(list.rhat, aes(value, fill=var)) +
            theme_linedraw() + 
            geom_density(alpha=0.5) +
            facet_wrap(~species)
g.rhat

list.rhat[list.rhat$value>1.2 & !is.na(list.rhat$value) & !is.na(list.rhat$value) &
            species=="Pink",]

# Plot: neff ===========================================================
g.neff <- ggplot(list.neff, aes(value, fill=var)) +
            theme_linedraw() + 
            geom_density(alpha=0.5) +
            facet_wrap(~species)
g.neff

# Plot: autocorr ===========================================================
g.ar <- ggplot(list.phi, aes(value, fill=var)) +
  theme_linedraw() + 
  geom_density(alpha=0.5) +
  facet_wrap(~species)
g.ar

# Plot: ratio ===========================================================

g.lims <- list.mu.ratio %>% group_by(region, species) %>% summarize(upper=quantile(value, 0.99))

g.ratio <- ggplot(list.mu.ratio, aes(x=region, y=exp(value), fill=var)) +
             scale_fill_colorblind() +
             theme_linedraw() +
             geom_violin() + 
             facet_wrap(~species) +
             # coord_flip(ylim=c(0,max(g.lims$upper)))
             coord_flip(ylim=c(0,5))

g.ratio












