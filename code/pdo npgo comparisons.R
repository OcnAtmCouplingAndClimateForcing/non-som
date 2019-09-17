# compare the distribution / SD / AR1 of 
# winter pdo and npgo between eras!

library(tidyverse)
library(zoo)

# download PDO / NPGO and process
download.file("http://jisao.washington.edu/pdo/PDO.latest", "~pdo")
names <- read.table("~pdo", skip=30, nrows=1, as.is = T)
pdo <- read.table("~pdo", skip=32, nrows=119, fill=T, col.names = names)
pdo$YEAR <- 1900:(1899+nrow(pdo)) # drop asterisks!
pdo <- pdo %>%
  gather(month, value, -YEAR) %>%
  arrange(YEAR)

download.file("http://www.oces.us/npgo/npgo.php", "~npgo")

npgo <- read.table("~npgo", skip=10, nrows=828, fill=T, col.names = c("Year", "month", "value"))

# calculate NDJFM means for each index
pdo$win.yr <- ifelse(pdo$month %in% c("NOV", "DEC"), pdo$YEAR+1, pdo$YEAR)

# limit to winter months only
win.pdo <- pdo %>%
  filter(month %in% c("NOV", "DEC", "JAN", "FEB", "MAR"))

win.pdo <- tapply(win.pdo$value, win.pdo$win.yr, mean)

npgo$win.yr <- ifelse(npgo$month %in% 11:12, npgo$Year+1, npgo$Year)

# limit to winter months only
win.npgo <- npgo %>%
  filter(month %in% c(11,12,1:3))

win.npgo <- tapply(win.npgo$value, win.npgo$win.yr, mean)

# and smoothed (2yr) values of each
npgo2 <- rollapply(win.npgo, 2, mean, align="right", fill=NA)
names(npgo2) <- 1950:2019

pdo2 <- rollapply(win.pdo, 2, mean, align="right", fill=NA)
names(pdo2) <- 1900:2019

dat <- data.frame(year=1950:2012,
                   pdo=win.pdo[names(win.pdo) %in% 1950:2012],
                   pdo2=pdo2[names(pdo2) %in% 1950:2012],
                   npgo=win.npgo[names(win.npgo) %in% 1950:2012],
                   npgo2=npgo2[names(npgo2) %in% 1950:2012])
plot <- dat %>%
  gather(index, value, -year)

plot$era <- as.factor(ifelse(plot$year <= 1988,1,2))

ggplot(plot, aes(value, fill=era)) +
  theme_linedraw() +
  geom_histogram(position="dodge", bins=12, color="black")+
  facet_wrap(~index)

ggplot(plot, aes(x=index, y=value, fill=era)) +
  theme_linedraw() +
  geom_boxplot(color="black")

ggsave("plots/pdo npgo boxplot.png")

# some other general plots
ts.plot <- data.frame(year=pdo$YEAR, month=1:12, dec.yr=NA, pdo=pdo$value, npgo=NA)
ts.plot$dec.yr <- ts.plot$year + (ts.plot$month-0.5)/12
ts.plot$npgo[601:1428] <- npgo$value

head(ts.plot)
ts.plot <- ts.plot %>%
  gather(index, value, -year, -month, -dec.yr)

ts.plot$color <- ifelse(ts.plot$value > 0, "red", "blue")
ts.plot$index <- reorder(ts.plot$index, desc(ts.plot$index))
ggplot(ts.plot, aes(dec.yr, value, fill=color)) +
  theme_linedraw() +
  geom_bar(stat="identity", color="black", size=0.05) +
  facet_wrap(~index, nrow=2, scales="free_y") +
  scale_fill_manual(values=c("blue", "red")) +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  xlim(1950,2018)

ggsave("plots/pdo npgo monthly time series.png")

npgo1 <- ts.plot %>%
  filter(index=="npgo", year %in% 1950:1988)
ar.npgo1 <- acf(npgo1$value)
plot(ar.npgo1)

npgo2 <- ts.plot %>%
  filter(index=="npgo", year %in% 1989:2012)
ar.npgo2 <- acf(npgo2$value)
ar.npgo2
plot(ar.npgo2)

pdo1 <- ts.plot %>%
  filter(index=="pdo", year %in% 1950:1988)
ar.pdo1 <- acf(pdo1$value)
ar.pdo1
plot(ar.pdo1)

pdo2 <- ts.plot %>%
  filter(index=="pdo", year %in% 1989:2012)
ar.pdo2 <- acf(pdo2$value)
ar.pdo2
plot(ar.pdo2)

# combine and plot
ar.plot <- data.frame(lag=rep(1:24, 4), ar=(c(ar.npgo1$acf[2:25], ar.npgo2$acf[2:25], ar.pdo1$acf[2:25], ar.pdo2$acf[2:25])),
                      index=rep(c("npgo", "pdo"), each=48), era=rep(c("1950-1988", "1989-2012", "1950-1988", "1989-2012"), each=24))

ar.plot$index <- reorder(ar.plot$index, desc(ar.plot$index))
ggplot(ar.plot, aes(lag, ar, fill=era)) +
  theme_linedraw() +
  geom_bar(stat="identity", position="dodge", color="black", size=0.05) +
  facet_wrap(~index, nrow=2) +
  ylab("autocorrelation") +
  theme(legend.position = c(0.8,0.9), legend.title = element_blank())

ggsave("plots/pdo npgo monthly ar by era.png")

# and finally, winters...

plot.win <- data.frame(year=1951:2018, pdo=win.pdo[names(win.pdo) %in% 1951:2018], npgo=win.npgo[names(win.npgo) %in% 1951:2018])

win.plot <- plot.win %>%
  gather(index, value, -year)

win.plot$color <- as.factor(ifelse(win.plot$value < 0, 1, 2))
win.plot$index <- reorder(win.plot$index, desc(win.plot$index))

ggplot(win.plot, aes(year, value, fill=color)) +
  theme_linedraw() +
  geom_bar(stat="identity", color="black", size=0.05) +
  facet_wrap(~index, scales="free_y", nrow=2) + 
  scale_fill_manual(values=c("blue", "red"), labels=c("1950-1988", "1989-2012")) +
  theme(legend.position = "none", axis.title.x = element_blank()) 

ggsave("plots/winter pdo npgo time series.png")

# finally, winter acf by era

npgo1 <- win.plot %>%
  filter(index=="npgo", year %in% 1950:1988)
ar.npgo1 <- acf(npgo1$value)
plot(ar.npgo1)

npgo2 <- win.plot %>%
  filter(index=="npgo", year %in% 1989:2012)
ar.npgo2 <- acf(npgo2$value)
ar.npgo2
plot(ar.npgo2)

pdo1 <- win.plot %>%
  filter(index=="pdo", year %in% 1950:1988)
ar.pdo1 <- acf(pdo1$value)
ar.pdo1
plot(ar.pdo1)

pdo2 <- win.plot %>%
  filter(index=="pdo", year %in% 1989:2012)
ar.pdo2 <- acf(pdo2$value)
ar.pdo2
plot(ar.pdo2)

# combine and plot
ar.plot <- data.frame(lag=rep(1:10, 4), ar=(c(ar.npgo1$acf[2:11], ar.npgo2$acf[2:11], ar.pdo1$acf[2:11], ar.pdo2$acf[2:11])),
                      index=rep(c("npgo", "pdo"), each=20), era=rep(c("1950-1988", "1989-2012", "1950-1988", "1989-2012"), each=10))

ar.plot$index <- reorder(ar.plot$index, desc(ar.plot$index))

ggplot(ar.plot, aes(as.factor(lag), ar, fill=era)) +
  theme_linedraw() +
  geom_bar(stat="identity", position="dodge", color="black", size=0.05) +
  facet_wrap(~index, nrow=2) +
  ylab("autocorrelation") +
  theme(legend.position = c(0.8,0.9), legend.title = element_blank()) +
  xlab("lag (years)")

ggsave("plots/pdo npgo winter ar by era.png")
