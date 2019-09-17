
library(MuMIn)
library(dplyr)
library(ggplot2)
library(nlme)
library(lmtest)
library(tidyverse)
library(ggpubr)
library(grid)

options(max.print = 99999)

#run.dat <- read.csv("/Users/MikeLitzow 1/Documents/R/coastwide-salmon/coastwide salmon data with different sst groupings.csv", row.names = 1)

#write.csv(run.dat, "salmon run dat.csv")

run.dat <- read.csv("data/salmon run dat.csv", row.names = 1)

# load pdo and npgo
download.file("http://jisao.washington.edu/pdo/PDO.latest", "~pdo")
names <- read.table("~pdo", skip=30, nrows=1, as.is = T)
pdo <- read.table("~pdo", skip=31, nrows=119, fill=T, col.names = names)
pdo$YEAR <- 1900:(1899+nrow(pdo)) # drop asterisks!
pdo <- pdo %>%
gather(month, value, -YEAR) %>%
arrange(YEAR)

# load npgo
download.file("http://www.oces.us/npgo/npgo.php", "~npgo")
npgo <- read.table("~npgo", skip=10, nrows=828, fill=T, col.names = c("Year", "month", "value"))

win.pdo <- pdo[pdo$month %in% c("NOV", "DEC", "JAN", "FEB", "MAR"),]
win.pdo$win.yr <- ifelse(win.pdo$month %in% c("NOV", "DEC"), win.pdo$YEAR+1, win.pdo$YEAR)
win.pdo <- tapply(win.pdo$value, win.pdo$win.yr, mean)

win.npgo <- npgo[npgo$month %in% c(11,12,1:3),]
win.npgo$win.yr <- ifelse(win.npgo$month %in% 11:12, win.npgo$Year+1, win.npgo$Year)

win.npgo <- tapply(win.npgo$value, win.npgo$win.yr, mean)

# smooth over various scales
# winter before and winter of ocean entry, winter of and winter after ocean entry,
# and winter before, winter of, and winter after ocean entry
npgo2 <- rollmean(win.npgo, 2, align="right", fill=NA)
npgo2a <- rollmean(win.npgo, 2, align="left", fill=NA)
npgo3 <- rollmean(win.npgo, 3, fill=NA)

pdo2 <- rollmean(win.pdo, 2, align="right", fill=NA)
pdo2a <- rollmean(win.pdo, 2, align="left", fill=NA)
pdo3 <- rollmean(win.pdo, 3, fill=NA)

run.dat$era <- as.factor(run.dat$era)

run.dat$pdo <- win.pdo[match(run.dat$entry.yr, names(win.pdo))]
run.dat$pdo2 <- pdo2[match(run.dat$entry.yr, names(pdo2))]
run.dat$pdo2a <- pdo2a[match(run.dat$entry.yr, names(pdo2a))]
run.dat$pdo3 <- pdo3[match(run.dat$entry.yr, names(pdo3))]

run.dat$npgo <- win.npgo[match(run.dat$entry.yr, names(win.npgo))]
run.dat$npgo2 <- npgo2[match(run.dat$entry.yr, names(npgo2))]
run.dat$npgo2a <- npgo2a[match(run.dat$entry.yr, names(npgo2a))]
run.dat$npgo3 <- npgo3[match(run.dat$entry.yr, names(npgo3))]

# now fit individual models
# using PDO/NPGO2a throughout (winter of and winter after ocean entry)

# first Pinks
form <-  as.formula(log(recruits/spawners) ~ 1 + spawners + pdo2a:era)

temp <- run.dat %>%
  filter(species=="Pink")

reps <- unique(temp$stock)

pink.early <- pink.late <- data.frame(stock=reps, est=NA, LCI=NA, UCI=NA)

for(i in 1:length(reps)){
  
  # fit to the pre-88/89 data
  mod <- gls(form, method = "REML", correlation=corAR1(), data = filter(temp, stock==reps[i]))
  
  pink.early$stock[i] <- reps[i] # re-define stock to make suire it's lining up correctly!
  
  pink.early$est[i] <- intervals(mod)$coef[3,2]
  pink.early$LCI[i] <- intervals(mod)$coef[3,1]
  pink.early$UCI[i] <- intervals(mod)$coef[3,3]
  
  # and the late era
  pink.late$stock[i] <- reps[i]
  
  pink.late$est[i] <- intervals(mod)$coef[4,2]
  pink.late$LCI[i] <- intervals(mod)$coef[4,1]
  pink.late$UCI[i] <- intervals(mod)$coef[4,3]
  
}

# combine results for each era
pink.early$era <- "Before 1988/89"
pink.late$era <- "After 1988/89"
pink.out <- rbind(pink.early, pink.late)

# add region
pink.out$region <- NA
for(i in 1:nrow(pink.out)){
  # i <- 1
  tt <- temp %>%
    filter(stock==pink.out$stock[i])
  pink.out$region[i] <- as.character(tt$region[1])
}


dodge <- position_dodge(width=0.9)

pink.out$stock <- reorder(pink.out$stock, rep(-pink.early$est, 2))
pink.out$plot.era <- as.factor(ifelse(pink.out$era=="Before 1988/89", 1, 2))

pink.out$region <- ifelse(pink.out$region=="EBS", "Eastern Bering Sea", pink.out$region)
pink.out$region <- ifelse(pink.out$region=="GOA", "Gulf of Alaska", pink.out$region)

# set pallette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
q <- ggplot(pink.out, aes(stock, est, fill=plot.era)) + 
  theme_linedraw() +
  geom_bar(position=dodge, stat="identity") + geom_hline(yintercept = 0, color="black", lwd=0.3) + 
  geom_errorbar(aes(ymax=UCI, ymin=LCI), position=dodge, width=0.3, size=0.2) +
  xlab("") + ylab("SST coefficient") + 
  facet_wrap(~region, scales="free") + 
  theme(legend.position = c(0.4,0.85), legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust=1)) +
  theme(legend.key.size =  unit(0.15, "in"), legend.text = element_text(size=8)) +
  scale_fill_manual(values=cb[c(2,6)], labels=c("Before 1988/89", "After 1988/89")) + 
  ggtitle("Pink - PDO")

# need to adjust widths of each panel
# Get the ggplot grob
gt = ggplotGrob(q)

# Check for the widths - you need to change the two that are set to 1null
gt$widths

# The required widths are 4 and 8

# Replace the default widths with relative widths:
gt$widths[5] = unit(0.5, "null")
gt$widths[9] = unit(1.5, "null")
gt$widths[13] = unit(0.8, "null")

# Draw the plot
png("plots/individual pink models with era-pdo interaction.png", 8, 3.8, units="in", res=300)
grid.newpage()
grid.draw(gt)
dev.off()

# pinks - npgo
form <-  as.formula(log(recruits/spawners) ~ 1 + spawners + npgo2a:era)

temp <- run.dat %>%
  filter(species=="Pink")

reps <- unique(temp$stock)

pink.early <- pink.late <- data.frame(stock=reps, est=NA, LCI=NA, UCI=NA)

for(i in 1:length(reps)){
  
  # fit to the pre-88/89 data
  mod <- gls(form, method = "REML", correlation=corAR1(), data = filter(temp, stock==reps[i]))
  
  pink.early$stock[i] <- reps[i] # re-define stock to make suire it's lining up correctly!
  
  pink.early$est[i] <- intervals(mod)$coef[3,2]
  pink.early$LCI[i] <- intervals(mod)$coef[3,1]
  pink.early$UCI[i] <- intervals(mod)$coef[3,3]
  
  # and the late era
  pink.late$stock[i] <- reps[i]
  
  pink.late$est[i] <- intervals(mod)$coef[4,2]
  pink.late$LCI[i] <- intervals(mod)$coef[4,1]
  pink.late$UCI[i] <- intervals(mod)$coef[4,3]
  
}

# combine results for each era
pink.early$era <- "Before 1988/89"
pink.late$era <- "After 1988/89"
pink.out <- rbind(pink.early, pink.late)

# add region
pink.out$region <- NA
for(i in 1:nrow(pink.out)){
  # i <- 1
  tt <- temp %>%
    filter(stock==pink.out$stock[i])
  pink.out$region[i] <- as.character(tt$region[1])
}


dodge <- position_dodge(width=0.9)

pink.out$stock <- reorder(pink.out$stock, rep(-pink.early$est, 2))
pink.out$plot.era <- as.factor(ifelse(pink.out$era=="Before 1988/89", 1, 2))

pink.out$region <- ifelse(pink.out$region=="EBS", "Eastern Bering Sea", pink.out$region)
pink.out$region <- ifelse(pink.out$region=="GOA", "Gulf of Alaska", pink.out$region)

# set pallette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
q <- ggplot(pink.out, aes(stock, est, fill=plot.era)) + 
  theme_linedraw() +
  geom_bar(position=dodge, stat="identity") + geom_hline(yintercept = 0, color="black", lwd=0.3) + 
  geom_errorbar(aes(ymax=UCI, ymin=LCI), position=dodge, width=0.3, size=0.2) +
  xlab("") + ylab("SST coefficient") + 
  facet_wrap(~region, scales="free") + 
  theme(legend.position = c(0.4,0.85), legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust=1)) +
  theme(legend.key.size =  unit(0.15, "in"), legend.text = element_text(size=8)) +
  scale_fill_manual(values=cb[c(2,6)], labels=c("Before 1988/89", "After 1988/89")) + 
  ggtitle("Pink - NPGO")

# need to adjust widths of each panel
# Get the ggplot grob
gt = ggplotGrob(q)

# Check for the widths - you need to change the two that are set to 1null
gt$widths

# The required widths are 4 and 8

# Replace the default widths with relative widths:
gt$widths[5] = unit(0.5, "null")
gt$widths[9] = unit(1.5, "null")
gt$widths[13] = unit(0.8, "null")

# Draw the plot
png("plots/individual pink models with era-npgo interaction.png", 8, 3.8, units="in", res=300)
grid.newpage()
grid.draw(gt)
dev.off()

################
# now sockeye
form <-  as.formula(log(recruits/spawners) ~ 1 + spawners + pdo2a:era)
run.dat$era <- as.factor(run.dat$era)

temp <- run.dat %>%
  filter(species=="Sockeye")

reps <- unique(temp$stock)

sock.early <- sock.late <- data.frame(stock=reps, est=NA, LCI=NA, UCI=NA)

for(i in 1:length(reps)){
  
  mod <- gls(form, method = "REML", correlation=corAR1(), data = filter(temp, stock==reps[i]))
  
  sock.early$stock[i] <- reps[i] # re-define stock to make suire it's lining up correctly!
  
  sock.early$est[i] <- intervals(mod)$coef[3,2]
  sock.early$LCI[i] <- intervals(mod)$coef[3,1]
  sock.early$UCI[i] <- intervals(mod)$coef[3,3]
  
  # and the late era
  sock.late$stock[i] <- reps[i]
  
  sock.late$est[i] <- intervals(mod)$coef[4,2]
  sock.late$LCI[i] <- intervals(mod)$coef[4,1]
  sock.late$UCI[i] <- intervals(mod)$coef[4,3]
  
}

# combine
sock.early$era <- "Before 1988/89"
sock.late$era <- "After 1988/89"
sock.out <- rbind(sock.early, sock.late)

# add region
sock.out$region <- NA
for(i in 1:nrow(sock.out)){
  # i <- 1
  tt <- temp %>%
    filter(stock==sock.out$stock[i])
  sock.out$region[i] <- as.character(tt$region[1])
}


dodge <- position_dodge(width=0.9)

sock.out$stock <- reorder(sock.out$stock, rep(-sock.early$est, 2))
sock.out$plot.era <- as.factor(ifelse(sock.out$era=="Before 1988/89", 1, 2))

sock.out$region <- ifelse(sock.out$region=="EBS", "Eastern Bering Sea", sock.out$region)
sock.out$region <- ifelse(sock.out$region=="GOA", "Gulf of Alaska", sock.out$region)

q <- ggplot(sock.out, aes(stock, est, fill=plot.era)) + 
  theme_linedraw() +
  geom_bar(position=dodge, stat="identity") + geom_hline(yintercept = 0, color="black", lwd=0.3) + 
  geom_errorbar(aes(ymax=UCI, ymin=LCI), position=dodge, width=0.3, size=0.2) +
  xlab("") + ylab("SST coefficient") + 
  facet_wrap(~region, scales="free") + 
  theme(legend.position = c(0.08, 0.13), legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust=1)) +
  theme(legend.key.size =  unit(0.15, "in"), legend.text = element_text(size=8)) +
  scale_fill_manual(values=cb[c(2,6)], labels=c("Before 1988/89", "After 1988/89")) + 
  ggtitle("Sockeye - PDO")


# Get the ggplot grob
gt = ggplotGrob(q)

# Check for the widths - you need to change the two that are set to 1null
gt$widths
# The required widths are 4 and 8

# Replace the default widths with relative widths:
gt$widths[5] = unit(0.7, "null")
gt$widths[9] = unit(1, "null")
gt$widths[13] = unit(0.7, "null")

# Draw the plot
png("plots/individual sockeye models with era-pdo interaction.png", 8, 3.8, units="in", res=300)
grid.newpage()
grid.draw(gt)
dev.off()

##
# sockeye - npgo
form <-  as.formula(log(recruits/spawners) ~ 1 + spawners + npgo2a:era)
run.dat$era <- as.factor(run.dat$era)

temp <- run.dat %>%
  filter(species=="Sockeye")

reps <- unique(temp$stock)

sock.early <- sock.late <- data.frame(stock=reps, est=NA, LCI=NA, UCI=NA)

for(i in 1:length(reps)){
  
  mod <- gls(form, method = "REML", correlation=corAR1(), data = filter(temp, stock==reps[i]))
  
  sock.early$stock[i] <- reps[i] # re-define stock to make suire it's lining up correctly!
  
  sock.early$est[i] <- intervals(mod)$coef[3,2]
  sock.early$LCI[i] <- intervals(mod)$coef[3,1]
  sock.early$UCI[i] <- intervals(mod)$coef[3,3]
  
  # and the late era
  sock.late$stock[i] <- reps[i]
  
  sock.late$est[i] <- intervals(mod)$coef[4,2]
  sock.late$LCI[i] <- intervals(mod)$coef[4,1]
  sock.late$UCI[i] <- intervals(mod)$coef[4,3]
  
}

# combine
sock.early$era <- "Before 1988/89"
sock.late$era <- "After 1988/89"
sock.out <- rbind(sock.early, sock.late)

# add region
sock.out$region <- NA
for(i in 1:nrow(sock.out)){
  # i <- 1
  tt <- temp %>%
    filter(stock==sock.out$stock[i])
  sock.out$region[i] <- as.character(tt$region[1])
}


dodge <- position_dodge(width=0.9)

sock.out$stock <- reorder(sock.out$stock, rep(-sock.early$est, 2))
sock.out$plot.era <- as.factor(ifelse(sock.out$era=="Before 1988/89", 1, 2))

sock.out$region <- ifelse(sock.out$region=="EBS", "Eastern Bering Sea", sock.out$region)
sock.out$region <- ifelse(sock.out$region=="GOA", "Gulf of Alaska", sock.out$region)

q <- ggplot(sock.out, aes(stock, est, fill=plot.era)) + 
  theme_linedraw() +
  geom_bar(position=dodge, stat="identity") + geom_hline(yintercept = 0, color="black", lwd=0.3) + 
  geom_errorbar(aes(ymax=UCI, ymin=LCI), position=dodge, width=0.3, size=0.2) +
  xlab("") + ylab("SST coefficient") + 
  facet_wrap(~region, scales="free") + 
  theme(legend.position = c(0.08, 0.13), legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust=1)) +
  theme(legend.key.size =  unit(0.15, "in"), legend.text = element_text(size=8)) +
  scale_fill_manual(values=cb[c(2,6)], labels=c("Before 1988/89", "After 1988/89")) + 
  ggtitle("Sockeye - NPGO")


# Get the ggplot grob
gt = ggplotGrob(q)

# Check for the widths - you need to change the two that are set to 1null
gt$widths
# The required widths are 4 and 8

# Replace the default widths with relative widths:
gt$widths[5] = unit(0.7, "null")
gt$widths[9] = unit(1, "null")
gt$widths[13] = unit(0.7, "null")

# Draw the plot
png("plots/individual sockeye models with era-npgo interaction.png", 8, 3.8, units="in", res=300)
grid.newpage()
grid.draw(gt)
dev.off()


################
# and chum
form <-  as.formula(log(recruits/spawners) ~ 1 + spawners + pdo2a:era)
run.dat$era <- as.factor(run.dat$era)

temp <- run.dat %>%
  filter(species=="Chum")

reps <- unique(temp$stock)

chum.early <- chum.late <- data.frame(stock=reps, est=NA, LCI=NA, UCI=NA)

for(i in 1:length(reps)){
  # i <- 1
  
  # fit to the pre-88/89 data
  mod <- gls(form, method = "REML", correlation=corAR1(), data = filter(temp, stock==reps[i]))
  
  chum.early$stock[i] <- reps[i] # re-define stock to make suire it's lining up correctly!
  
  chum.early$est[i] <- intervals(mod)$coef[3,2]
  chum.early$LCI[i] <- intervals(mod)$coef[3,1]
  chum.early$UCI[i] <- intervals(mod)$coef[3,3]
  
  # and the late era
  chum.late$stock[i] <- reps[i]
  
  chum.late$est[i] <- intervals(mod)$coef[4,2]
  chum.late$LCI[i] <- intervals(mod)$coef[4,1]
  chum.late$UCI[i] <- intervals(mod)$coef[4,3]
  
}

# combine
chum.early$era <- "Before 1988/89"
chum.late$era <- "After 1988/89"
chum.out <- rbind(chum.early, chum.late)

# add region
chum.out$region <- NA
for(i in 1:nrow(chum.out)){
  # i <- 1
  tt <- temp %>%
    filter(stock==chum.out$stock[i])
  chum.out$region[i] <- as.character(tt$region[1])
}


dodge <- position_dodge(width=0.9)

chum.out$stock <- reorder(chum.out$stock, rep(-chum.early$est, 2))
chum.out$plot.era <- as.factor(ifelse(chum.out$era=="Before 1988/89", 1, 2))

chum.out$region <- ifelse(chum.out$region=="EBS", "Eastern Bering Sea", chum.out$region)
chum.out$region <- ifelse(chum.out$region=="GOA", "Gulf of Alaska", chum.out$region)

lm <- 0.03

q <- ggplot(chum.out, aes(stock, est, fill=plot.era)) + 
  theme_linedraw() +
  geom_bar(position=dodge, stat="identity") + geom_hline(yintercept = 0, color="black", lwd=0.3) + 
  geom_errorbar(aes(ymax=UCI, ymin=LCI), position=dodge, width=0.3, size=0.2) +
  xlab("") + ylab("SST coefficient") + 
  facet_wrap(~region, scales="free") + 
  theme(legend.position = c(0.5, 0.09), legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust=1)) +
  theme(legend.key.size =  unit(0.15, "in"), legend.text = element_text(size=8), legend.margin = margin(lm,lm,lm,lm,"cm")) +
  scale_fill_manual(values=cb[c(2,6)], labels=c("Before 1988/89", "After 1988/89")) + 
  ggtitle("Chum - PDO")


# Get the ggplot grob
gt = ggplotGrob(q)

# Check for the widths - you need to change the two that are set to 1null
gt$widths
# The required widths are 4 and 8

# Replace the default widths with relative widths:
gt$widths[5] = unit(0.4, "null")
gt$widths[9] = unit(1.1, "null")
gt$widths[13] = unit(0.8, "null")

# Draw the plot
png("plots/individual chum models with era-pdo interaction.png", 8, 3.8, units="in", res=300)
grid.newpage()
grid.draw(gt)
dev.off()

# chum - npgo
form <-  as.formula(log(recruits/spawners) ~ 1 + spawners + npgo2a:era)
run.dat$era <- as.factor(run.dat$era)

temp <- run.dat %>%
  filter(species=="Chum")

reps <- unique(temp$stock)

chum.early <- chum.late <- data.frame(stock=reps, est=NA, LCI=NA, UCI=NA)

for(i in 1:length(reps)){
  # i <- 1
  
  # fit to the pre-88/89 data
  mod <- gls(form, method = "REML", correlation=corAR1(), data = filter(temp, stock==reps[i]))
  
  chum.early$stock[i] <- reps[i] # re-define stock to make suire it's lining up correctly!
  
  chum.early$est[i] <- intervals(mod)$coef[3,2]
  chum.early$LCI[i] <- intervals(mod)$coef[3,1]
  chum.early$UCI[i] <- intervals(mod)$coef[3,3]
  
  # and the late era
  chum.late$stock[i] <- reps[i]
  
  chum.late$est[i] <- intervals(mod)$coef[4,2]
  chum.late$LCI[i] <- intervals(mod)$coef[4,1]
  chum.late$UCI[i] <- intervals(mod)$coef[4,3]
  
}

# combine
chum.early$era <- "Before 1988/89"
chum.late$era <- "After 1988/89"
chum.out <- rbind(chum.early, chum.late)

# add region
chum.out$region <- NA
for(i in 1:nrow(chum.out)){
  # i <- 1
  tt <- temp %>%
    filter(stock==chum.out$stock[i])
  chum.out$region[i] <- as.character(tt$region[1])
}


dodge <- position_dodge(width=0.9)

chum.out$stock <- reorder(chum.out$stock, rep(-chum.early$est, 2))
chum.out$plot.era <- as.factor(ifelse(chum.out$era=="Before 1988/89", 1, 2))

chum.out$region <- ifelse(chum.out$region=="EBS", "Eastern Bering Sea", chum.out$region)
chum.out$region <- ifelse(chum.out$region=="GOA", "Gulf of Alaska", chum.out$region)

lm <- 0.03

q <- ggplot(chum.out, aes(stock, est, fill=plot.era)) + 
  theme_linedraw() +
  geom_bar(position=dodge, stat="identity") + geom_hline(yintercept = 0, color="black", lwd=0.3) + 
  geom_errorbar(aes(ymax=UCI, ymin=LCI), position=dodge, width=0.3, size=0.2) +
  xlab("") + ylab("SST coefficient") + 
  facet_wrap(~region, scales="free") + 
  theme(legend.position = c(0.5, 0.09), legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust=1)) +
  theme(legend.key.size =  unit(0.15, "in"), legend.text = element_text(size=8), legend.margin = margin(lm,lm,lm,lm,"cm")) +
  scale_fill_manual(values=cb[c(2,6)], labels=c("Before 1988/89", "After 1988/89")) + 
  ggtitle("Chum - NPGO")


# Get the ggplot grob
gt = ggplotGrob(q)

# Check for the widths - you need to change the two that are set to 1null
gt$widths
# The required widths are 4 and 8

# Replace the default widths with relative widths:
gt$widths[5] = unit(0.4, "null")
gt$widths[9] = unit(1.1, "null")
gt$widths[13] = unit(0.8, "null")

# Draw the plot
png("plots/individual chum models with era-npgo interaction.png", 8, 3.8, units="in", res=300)
grid.newpage()
grid.draw(gt)
dev.off()