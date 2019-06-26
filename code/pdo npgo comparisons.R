# compare the distribution / SD / AR1 of 
# winter pdo and npgo between eras!

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
pdo <- pdo %>%
  filter(month %in% c("NOV", "DEC", "JAN", "FEB", "MAR"))

win.pdo <- tapply(pdo$value, pdo$win.yr, mean)

npgo$win.yr <- ifelse(npgo$month %in% 11:12, npgo$Year+1, npgo$Year)

# limit to winter months only
npgo <- npgo %>%
  filter(month %in% c(11,12,1:3))

win.npgo <- tapply(npgo$value, npgo$win.yr, mean)

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

