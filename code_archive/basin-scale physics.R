library(ncdf4)
library(zoo)
library(gplots)
library(dplyr)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(tidyr)
library(nlme)
library(pracma)
library(ggplot2)
library(car)
library(FactoMineR)
library(ggpubr)
library(mgcv)


# using monthly NCEP/NCAR!
# these data are in the google drive data folder
# nc.slp <- nc_open("/Users/MikeLitzow 1/Documents/R/climate scripts/ncep_ncar_monthly_slp.nc")
# 
# # now process SLP data - first, extract dates
# raw <- ncvar_get(nc.slp, "TIME")  # seconds since 1-1-1970
# # h <- raw/(24*60*60)
# d <- dates(raw, origin = c(1,1,0001))
# yr <- years(d)
# 
# x <- ncvar_get(nc.slp, "LON53_101")
# y <- ncvar_get(nc.slp, "LAT45_69")

nc.slp <- nc_open("/Users/MikeLitzow 1/Documents/R/climate scripts/monthly.SLP.full.N.Pac.NCEP.NCAR.nc")

# now process SLP data - first, extract dates
raw <- ncvar_get(nc.slp, "TIME")  # seconds since 1-1-1970
# h <- raw/(24*60*60)
d <- dates(raw, origin = c(1,1,0001))
m <- months(d)
yr <- years(d)
dec.yr <- as.numeric(as.character(yr)) + (as.numeric(m)-0.5)/12
# and lat/long
x <- ncvar_get(nc.slp, "LON53_101")
y <- ncvar_get(nc.slp, "LAT45_65")
SLP <- ncvar_get(nc.slp, "SLP", verbose = F)
# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SLP <- aperm(SLP, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SLP <- matrix(SLP, nrow=dim(SLP)[1], ncol=prod(dim(SLP)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   
dimnames(SLP) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# load pdo
download.file("http://jisao.washington.edu/pdo/PDO.latest", "~pdo")
names <- read.table("~pdo", skip=30, nrows=1, as.is = T)
pdo <- read.table("~pdo", skip=31, nrows=119, fill=T, col.names = names)
pdo$YEAR <- 1900:(1899+nrow(pdo)) # drop asterisks!
pdo <- pdo %>%
  gather(month, value, -YEAR) %>%
  arrange(YEAR)

# load npgo
# download.file("http://www.oces.us/npgo/npgo.php", "~npgo")
npgo <- read.table("~npgo", skip=10, nrows=828, fill=T, col.names = c("Year", "month", "value"))




###
# and SLP-PDO
# reload SLP to make things easy...


# # save to plot below
# x.slp <- x
# y.slp <- y

# SLP <- ncvar_get(nc.slp, "SLP", verbose = F)
# # Change data from a 3-D array to a matrix of monthly data by grid point:
# # First, reverse order of dimensions ("transpose" array)
# SLP <- aperm(SLP, 3:1)  
# 
# # Change to matrix with column for each grid point, rows for monthly means
# SLP <- matrix(SLP, nrow=dim(SLP)[1], ncol=prod(dim(SLP)[2:3]))  
# 
# # Keep track of corresponding latitudes and longitudes of each column:
# lat <- rep(y, length(x))   
# lon <- rep(x, each = length(y))   
# dimnames(SLP) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# # load pdo
# names <- read.table("~pdo", skip=30, nrows=1, as.is = T)
# pdo <- read.table("~pdo", skip=31, nrows=119, fill=T, col.names = names)
# pdo$YEAR <- 1900:(1899+nrow(pdo)) # drop asterisks!
# pdo <- pdo %>%
#   gather(month, value, -YEAR) %>%
#   arrange(YEAR)
# 
# # load npgo
# # download.file("http://www.oces.us/npgo/npgo.php", "~npgo")
# npgo <- read.table("~npgo", skip=10, nrows=828, fill=T, col.names = c("Year", "month", "value"))

# smooth SLP with three month rolling mean

# limit SLP to NDJ
SLP.NDJ <- SLP[m %in% c("Nov", "Dec", "Jan"),]

yr <- as.numeric(as.character(yr))
win.yr <- ifelse(m %in% c("Nov", "Dec"), yr+1, yr)
win.yr <- win.yr[m %in% c("Nov", "Dec", "Jan")]

# get NDJ means for each cell 
rownames(SLP.NDJ) # 1949-2019 are complete!

ff <- function(x) tapply(x, win.yr, mean)
SLP.NDJ <- apply(SLP.NDJ, 2, ff)

# limit pdo and npgo to FMA
pdo <- pdo %>%
  filter(month %in% c("FEB", "MAR", "APR"))

PDO.FMA <- tapply(pdo$value, pdo$YEAR, mean) # complete through 2018

npgo <- npgo %>%
  filter(month %in% 2:4)

NPGO.FMA <- tapply(npgo$value, npgo$Year, mean) # complete through 2018

# separate SLP, PDO, and NPGO into era-specific chunks for regression maps

SLP1 <- SLP.NDJ[rownames(SLP.NDJ) %in% 1950:1988,]
SLP2 <- SLP.NDJ[rownames(SLP.NDJ) %in% 1989:2012,]

PDO1 <- PDO.FMA[names(PDO.FMA) %in% 1950:1988]
PDO2 <- PDO.FMA[names(PDO.FMA) %in% 1989:2012]

NPGO1 <- NPGO.FMA[names(NPGO.FMA) %in% 1950:1988]
NPGO2 <- NPGO.FMA[names(NPGO.FMA) %in% 1989:2012]

# separate regressions in each era!
pdo.regr1 <- pdo.regr2 <- npgo.regr1 <- npgo.regr2 <- NA

for(i in 1:ncol(SLP1)){
  #  i <- 1
  mod <- lm(SLP1[,i] ~ PDO1)
  pdo.regr1[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SLP2[,i] ~ PDO2)
  pdo.regr2[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SLP1[,i] ~ NPGO1)
  npgo.regr1[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP2[,i] ~ NPGO2)
  npgo.regr2[i] <- summary(mod)$coef[2,1] 
}

# calculate era differences for each
pdo.diff <- pdo.regr2 - pdo.regr1
npgo.diff <- npgo.regr2 - npgo.regr1
diff.lim <- range(pdo.diff, npgo.diff)

# # and plot
# # set up color schemes
# new.col <- my.col <- tim.colors(64)
# grays <- c("gray98", "gray97", "gray96", "gray95", "gray94", "gray93", "gray92", "gray91", "gray90", "gray89", "gray88",
#            "gray87", "gray86", "gray85", "gray84", "gray83", "gray82", "gray81", "gray80", "gray79", "gray78", "gray77",
#            "gray76", "gray75", "gray74", "gray73", "gray72", "gray71", "gray71", "gray71", "gray70", "gray69", "gray68")
# 
# 
# my.col[1:33] <- grays
# my.col[22:43] <- c(grays[11:1], grays[1:11])
# 
# new.col[27:36] <- c(grays[5:1], grays[1:5])
# 
# lim <- range(pdo.regr1, pdo.regr2, npgo.regr1, npgo.regr2)
# 
# 
# # setup the layout
# mt.cex <- 1.1
# l.mar <- 3
# l.cex <- 0.8
# l.l <- 0.2
# tc.l <- -0.2


###########################
# get era regression coefficients for each square

# define the squares
xx1 <- c(186, 186, 201.5, 201.5, 186)
yy1 <- c(51.5, 56.5, 56.5, 51.5, 51.5)

xx2 <- c(219, 219, 231.5, 231.5, 219)
yy2 <- c(39, 49, 49, 39, 39)

xx3 <- c(183.5, 183.5, 204, 204, 183.5)
yy3 <- c(46.5, 56.5, 56.5, 46.5, 46.5)


PDOkeep1 <- PDOkeep2 <- NPGOkeep <- NA
pdoSLP1.era1 <- pdoSLP2.era1 <- npgoSLP.era1 <- SLP1
pdoSLP1.era2 <- pdoSLP2.era2 <- npgoSLP.era2 <- SLP2

for(i in 1:length(lat)){
  # i <- 1
  PDOkeep1[i] <- inpolygon(lon[i], lat[i], xx1, yy1)
  PDOkeep2[i] <- inpolygon(lon[i], lat[i], xx2, yy2)
  NPGOkeep[i] <- inpolygon(lon[i], lat[i], xx3, yy3)
}

pdoSLP1.era1[,!PDOkeep1] <- NA
pdoSLP1.era2[,!PDOkeep1] <- NA
pdoSLP2.era1[,!PDOkeep2] <- NA
pdoSLP2.era2[,!PDOkeep2] <- NA
npgoSLP.era1[,!NPGOkeep] <- NA
npgoSLP.era2[,!NPGOkeep] <- NA

pdoSLP1.era1 <- rowMeans(pdoSLP1.era1, na.rm=T)
pdoSLP1.era2 <- rowMeans(pdoSLP1.era2, na.rm=T)

pdoSLP2.era1 <- rowMeans(pdoSLP2.era1, na.rm=T)
pdoSLP2.era2 <- rowMeans(pdoSLP2.era2, na.rm=T)

npgoSLP.era1 <- rowMeans(npgoSLP.era1, na.rm=T)
npgoSLP.era2 <- rowMeans(npgoSLP.era2, na.rm=T)

plot(pdoSLP1.era1, PDO1, pch=19, col="red")
points(pdoSLP1.era2, PDO2, pch=19, col="blue")
abline(lsfit(pdoSLP1.era1, PDO1), col="red")
abline(lsfit(pdoSLP1.era2, PDO2), col="blue")

plot(pdoSLP2.era1, PDO1, pch=19, col="red")
points(pdoSLP2.era2, PDO2, pch=19, col="blue")
abline(lsfit(pdoSLP2.era1, PDO1), col="red")
abline(lsfit(pdoSLP2.era2, PDO2), col="blue")

plot(npgoSLP.era1, NPGO1, pch=19, col="red", xlim=range(npgoSLP.era1, npgoSLP.era2), ylim=range(NPGO1, NPGO2))
points(npgoSLP.era2, NPGO2, pch=19, col="blue")
abline(lsfit(npgoSLP.era1, NPGO1), col="red")
abline(lsfit(npgoSLP.era2, NPGO2), col="blue")
# 
# # make regression maps for each era and difference
# 
# # setup the layout
# mt.cex <- 1.1
# l.mar <- 3
# l.cex <- 0.8
# l.l <- 0.2
# tc.l <- -0.2
# 
# cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# new.col <- my.col <- tim.colors(64)
# grays <- c("gray98", "gray98", "gray97", "gray97", "gray96", "gray96", "gray95", "gray95", "gray94", "gray94", "gray93",
#            "gray93", "gray92", "gray92", "gray91", "gray91", "gray90", "gray90", "gray89", "gray89", "gray88", "gray88",
#            "gray87", "gray87", "gray86", "gray86", "gray85", "gray85", "gray84", "gray84", "gray83", "gray83", "gray82")
# 
# 
# my.col[1:33] <- grays
# # my.col[22:43] <- c(grays[11:1], grays)
# new.col[27:36] <- c(grays[5:1], grays[1:5])
# 
# # xlim <- c(min(dec.yr), max(dec.yr.t))
# 
# par(mar=c(1.25,1.25,1.25,1),  tcl=tc.l, mgp=c(1.5,0.3,0), las=1, mfrow=c(4,4), cex.axis=0.8, cex.lab=0.8, oma=c(0,2,2,0))
# # now atmospheric forcing of PDO/NPGO
# lim <- range(pdo.regr1, pdo.regr2, npgo.regr1, npgo.regr2)
# diff.lim <- range(pdo.diff, npgo.diff)
# # PDO first era
# z <- pdo.regr1  # replace elements NOT corresponding to land with loadings!
# z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
# #image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
# #           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
# image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#            xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
# 
# contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
# map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
# map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
# map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
# map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
# map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
# map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
# map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
# mtext("d", adj=0.05, line=-1.4, cex=1)
# 
# xx1 <- c(183.5, 183.5, 204, 204, 183.5)
# yy1 <- c(46.5, 54, 54, 46.5, 46.5)
# 
# xx2 <- c(196, 196, 216.5, 216.5, 196)
# yy2 <- c(39, 46.5, 46.5, 39, 39)
# 
# lines(xx1, yy1, lwd=1.5, col="magenta")
# lines(xx2, yy2, lwd=1.5, col="magenta")
# 
# # PDO second era
# z <- pdo.regr2  # replace elements NOT corresponding to land with loadings!
# z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
# #image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
# #           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
# image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#            xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
# 
# contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
# map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
# map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
# map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
# map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
# map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
# map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
# map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
# mtext("e", adj=0.05, line=-1.4, cex=1)
# 
# lines(xx1, yy1, lwd=1.5, col="magenta")
# lines(xx2, yy2, lwd=1.5, col="magenta")
# 
# par(mar=c(1.5,3,1.5,1))
# 
# plot(dat$year, dat$pdoNS, type="l", xlab="", ylab="North:South forcing ratio", xlim=xlim, col=cb[6])
# abline(h=mean(dat$pdoNS, na.rm=T))
# abline(v=1989, lty=2)
# mtext("f", adj=0.05, line=-1.4, cex=1)
# 
# par(mar=c(1.25,1.25,1.25,1))
# 
# # npgo first era
# z <- npgo.regr1  # replace elements NOT corresponding to land with loadings!
# z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
# # image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
# #            xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
# image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#            xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
# 
# contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
# map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
# map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
# map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
# map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
# map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
# map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
# map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
# mtext("g", adj=0.05, line=-1.4, cex=1)
# 
# xx3 <- c(194, 194, 199, 211, 216.5, 216.5, 194)
# yy3 <- c(51.5, 54, 56.5, 61.5, 61.5, 51.5, 51.5)
# 
# lines(xx3, yy3, lwd=1.5, col="magenta")
# 
# # npgo second era
# z <- npgo.regr2  # replace elements NOT corresponding to land with loadings!
# z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
# # image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
# #            xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
# image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#            xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
# 
# contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
# map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
# map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
# map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
# map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
# map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
# map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
# map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
# mtext("h", adj=0.05, line=-1.4, cex=1)
# 
# dat$pdoNS <- dat$pdo.regr1/dat$pdo.regr2
#   
# plot.dat$SLP.PDO.NS  <- dat$pdoNS[match(plot.dat$dec.yr, dat$year)]
# plot.dat$SLP.NPGO  <- dat$npgo.regr[match(plot.dat$dec.yr, dat$year)]
# 
# plot.dat <- dat %>%
#   select(year, pdo.regr1, pdo.regr2, npgo.regr) %>%
#   gather(key, value, -year)
# 
# ggplot(plot.dat, aes(year, value, color=key)) +
#   theme_linedraw() +
#   geom_line() +
#   ylim(c(-0.5, -2.8)) + 
#   xlim(c(1962,2005))
# 
# dat$pdo.ratio <- ifelse(is.na(dat$pdo.regr1),NA, dat$pdo.regr1/dat$pdo.regr2)
# 
# plot.dat <- dat %>%
#   select(year, pdo.ratio, npgo.regr) %>%
#   gather(key, value, -year)
# 
# ggplot(plot.dat, aes(year, value)) +
#   theme_linedraw() +
#   geom_line() +
#   facet_wrap(~key, scales="free") + 
#   xlim(c(1962,2005))



#######
# combined plot
png("slp.pdo.npgo.maps.png", 8,6, units="in", res=300)

# setup the layout
mt.cex <- 1.1
l.mar <- 3
l.cex <- 0.8
l.l <- 0.2
tc.l <- -0.2

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
new.col <- my.col <- tim.colors(64)
grays <- c("gray98", "gray98", "gray97", "gray97", "gray96", "gray96", "gray95", "gray95", "gray94", "gray94", "gray93",
           "gray93", "gray92", "gray92", "gray91", "gray91", "gray90", "gray90", "gray89", "gray89", "gray88", "gray88",
           "gray87", "gray87", "gray86", "gray86", "gray85", "gray85", "gray84", "gray84", "gray83", "gray83", "gray82")


my.col[1:33] <- grays
# my.col[22:43] <- c(grays[11:1], grays)
new.col[27:36] <- c(grays[5:1], grays[1:5])

xlim <- c(160,250)

par(mar=c(1.25,1.25,1.25,1),  tcl=tc.l, mgp=c(1.5,0.3,0), las=1, mfrow=c(2,3), cex.axis=0.8, cex.lab=0.8, oma=c(0,2,2,0))



lim <- range(pdo.regr1, pdo.regr2, npgo.regr1, npgo.regr2)
# PDO first era
z <- pdo.regr1  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
#image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68), xlim=xlim,
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
#mtext("a", adj=0.05, line=-1.4, cex=1)
mtext("a) SLP-PDO 1950-1988", adj=0)


# PDO second era
z <- pdo.regr2  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
#image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68), xlim=xlim,
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
# mtext("b", adj=0.05, line=-1.4, cex=1)
mtext("b) SLP-PDO 1989-2012", adj=0)

# PDO diff
z <- pdo.diff  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
#image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x,y,z, col=new.col, zlim=c(-diff.lim[2], diff.lim[2]), ylim=c(20,68), xlim=xlim,
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
# mtext("c", adj=0.05, line=-1.4, cex=1)
mtext("c) SLP-PDO: difference", adj=0)

# xx1 <- c(186, 186, 201.5, 201.5, 186)
# yy1 <- c(51.5, 56.5, 56.5, 51.5, 51.5)
# 
# xx2 <- c(219, 219, 231.5, 231.5, 219)
# yy2 <- c(39, 49, 49, 39, 39)
# # 
# lines(xx1, yy1, lwd=1.5, col="magenta")
# lines(xx2, yy2, lwd=1.5, col="magenta")

# npgo first era
z <- npgo.regr1  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
# image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#            xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68), xlim=xlim,
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
# mtext("d", adj=0.05, line=-1.4, cex=1)
mtext("d) SLP-NPGO 1950-1988", adj=0)
# npgo second era
z <- npgo.regr2  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
# image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#            xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68), xlim=xlim,
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
# mtext("e", adj=0.05, line=-1.4, cex=1)
mtext("e) SLP-NPGO 1989-2012", adj=0)

# npgo diff
z <- npgo.diff # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
# image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#            xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x,y,z, col=new.col, zlim=c(-diff.lim[2], diff.lim[2]), ylim=c(20,68), xlim=xlim,
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
# mtext("f", adj=0.05, line=-1.4, cex=1)
mtext("f) SLP-NPGO: difference", adj=0)
# xx3 <- c(194, 194, 199, 211, 216.5, 216.5, 194)
# yy3 <- c(51.5, 54, 56.5, 61.5, 61.5, 51.5, 51.5)
# # 
# xx3 <- c(183.5, 183.5, 204, 204, 183.5)
# yy3 <- c(46.5, 56.5, 56.5, 46.5, 46.5)
# # 
# # 
# lines(xx3, yy3, lwd=1.5, col="magenta")

dev.off()

par(mar=c(1.5,3,1.5,1))

plot(dat$year, dat$npgo.regr, type="l", xlab="", ylab="Regression coef (Pa)", xlim=xlim, col=cb[6],
     ylim=c(max(dat$npgo.regr, na.rm=T), min(dat$npgo.regr, na.rm=T)))
abline(h=mean(dat$npgo.regr, na.rm=T))
abline(v=1989, lty=2)
mtext("i", adj=0.05, line=-1.4, cex=1)

par(mar=c(1.25,1.25,1.25,1))


# finally, the NPGO-SST plots
lim <- range(npgo.sst.regr1, npgo.sst.regr2)
z <- rep(NA, ncol(SST))
# z[!land] <- npgo.pattern1  # replace elements NOT corresponding to land with loadings!
# z[!land] <- npgo.eof.r1  # replace elements NOT corresponding to land with loadings!
z[!land] <- npgo.sst.regr1  # replace elements NOT corresponding to land with loadings!

z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=new.col, zlim=c(lim[1], -lim[1]),
           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))


contour(x.t,y.t,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("j", adj=0.05, line=-1.4, cex=1)

z <- rep(NA, ncol(SST))
# z[!land] <- npgo.pattern2  # replace elements NOT corresponding to land with loadings!
# z[!land] <- npgo.eof.r2
z[!land] <- npgo.sst.regr2
z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=new.col, zlim=c(lim[1], -lim[1]),
           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))


contour(x.t,y.t,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("k", adj=0.05, line=-1.4, cex=1)

par(mar=c(1.5,3,1.5,1))

plot(plot.dat$dec.yr, plot.dat$PDO.NPGO.cor, type="l", xlab="", ylab="PDO-NPGO correlation", xlim=xlim, col=cb[6])
abline(h=0)
abline(v=1989, lty=2)
mtext("l", adj=0.05, line=-1.4, cex=1)
dev.off()


########
# older version
png("older larger basin scale combined plot.png", 7,10, units="in", res=300)

# setup the layout
mt.cex <- 1.1
l.mar <- 3
l.cex <- 0.8
l.l <- 0.2
tc.l <- -0.2

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
 my.col <- tim.colors(64)
grays <- c("gray98", "gray98", "gray97", "gray97", "gray96", "gray96", "gray95", "gray95", "gray94", "gray94", "gray93",
           "gray93", "gray92", "gray92", "gray91", "gray91", "gray90", "gray90", "gray89", "gray89", "gray88", "gray88",
           "gray87", "gray87", "gray86", "gray86", "gray85", "gray85", "gray84", "gray84", "gray83", "gray83", "gray82")


my.col[1:33] <- grays
# my.col[22:43] <- c(grays[11:1], grays)
# new.col[27:36] <- c(grays[5:1], grays[1:5])

xlim <- c(min(dec.yr), max(dec.yr.t))

par(mar=c(1.5,1.5,1.5,1),  tcl=tc.l, mgp=c(1.5,0.3,0), las=1, mfrow=c(5,3), cex.axis=0.8, cex.lab=0.8, oma=c(0,2,2,0))

lim <- range(SSTanom1, SSTanom2, na.rm=T)

z <- SSTanom1   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=tim.colors(64), xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim[2], lim[2]), 
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x.t,y.t,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,80),add=T, lwd=1)
mtext("a", adj=0.05, line=-1.4, cex=1.1)

z <- SSTanom2   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=tim.colors(64), xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim[2], lim[2]),
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x.t,y.t,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,80),add=T, lwd=1)
mtext("b", adj=0.05, line=-1.4, cex=1.1)

par(mar=c(1.5,3,1.5,1))
plot(dec.yr.t, SST.anomTS, type="l", xlab="", ylab="ºC wrt 1951-1980", xlim=xlim, col=cb[6])
abline(h=0)
abline(v=1989, lty=2)
mtext("c", adj=0.05, line=-1.4, cex=1.1)

par(mar=c(1.5,1.5,1.5,1))

###
# now AL SD

lim <- range(SLPsd1, SLPsd2)

xx <- c(186.25, 186.25, 206, 206, 186.25)
yy <- c(46.25, 56.25, 56.25, 46.25, 46.25)


z <- SLPsd1   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image.plot(x,y,z, col=my.col, xlab = "", ylab = "", zlim=lim, ylim=c(20,70), yaxt="n", xaxt="n",
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,80),add=T, lwd=1)
lines(xx,yy, lwd=2, col="magenta")
mtext("d", adj=0.05, line=-1.4, cex=1.1)

z <- SLPsd2   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image.plot(x,y,z, col=my.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=lim, ylim=c(20,70), yaxt="n", xaxt="n", 
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,80),add=T, lwd=1)
lines(xx,yy, lwd=2, col="magenta")
mtext("e", adj=0.05, line=-1.4, cex=1.1)

par(mar=c(1.5,3,1.5,1))

plot(dec.yr, SLP.sd, type="l", xlab="", ylab="Standard dev. (Pa)", xlim=xlim, col=cb[6])
abline(h=mean(SLP.sd, na.rm=T))
abline(v=1989, lty=2)
mtext("f", adj=0.05, line=-1.4, cex=1.1)

par(mar=c(1.5,1.5,1.5,1))

# now atmospheric forcing of PDO/NPGO
lim <- range(pdo.regr1, pdo.regr2, npgo.regr1, npgo.regr2)
# PDO first era
z <- pdo.regr1  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.slp)))  # Convert vector to matrix and transpose for plotting
#image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x.slp,y.slp,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
      xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x.slp,y.slp,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("g", adj=0.05, line=-1.4, cex=1.1)

xx1 <- c(183.5, 183.5, 204, 204, 183.5)
yy1 <- c(46.5, 54, 54, 46.5, 46.5)

xx2 <- c(196, 196, 216.5, 216.5, 196)
yy2 <- c(39, 46.5, 46.5, 39, 39)

lines(xx1, yy1, lwd=1.5, col="magenta")
lines(xx2, yy2, lwd=1.5, col="magenta")

# PDO second era
z <- pdo.regr2  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.slp)))  # Convert vector to matrix and transpose for plotting
#image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x.slp,y.slp,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x.slp,y.slp,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("h", adj=0.05, line=-1.4, cex=1.1)

lines(xx1, yy1, lwd=1.5, col="magenta")
lines(xx2, yy2, lwd=1.5, col="magenta")

par(mar=c(1.5,3,1.5,1))

plot(dat$year, dat$pdoNS, type="l", xlab="", ylab="North:South forcing ratio", xlim=xlim, col=cb[6])
abline(h=mean(dat$pdoNS, na.rm=T))
abline(v=1989, lty=2)
mtext("i", adj=0.05, line=-1.4, cex=1.1)

par(mar=c(1.5,1.5,1.5,1))

# npgo first era
z <- npgo.regr1  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.slp)))  # Convert vector to matrix and transpose for plotting
# image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#            xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x.slp,y.slp,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x.slp,y.slp,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("j", adj=0.05, line=-1.4, cex=1.1)

xx3 <- c(194, 194, 199, 211, 216.5, 216.5, 194)
yy3 <- c(51.5, 54, 56.5, 61.5, 61.5, 51.5, 51.5)

lines(xx3, yy3, lwd=1.5, col="magenta")

# npgo second era
z <- npgo.regr2  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.slp)))  # Convert vector to matrix and transpose for plotting
# image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#            xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x.slp,y.slp,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x.slp,y.slp,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("k", adj=0.05, line=-1.4, cex=1.1)

xx3 <- c(194, 194, 199, 211, 216.5, 216.5, 194)
yy3 <- c(51.5, 54, 56.5, 61.5, 61.5, 51.5, 51.5)

lines(xx3, yy3, lwd=1.5, col="magenta")

par(mar=c(1.5,3,1.5,1))

plot(dat$year, dat$npgo.regr, type="l", xlab="", ylab="Regression coef (Pa)", xlim=xlim, col=cb[6],
     ylim=c(max(dat$npgo.regr, na.rm=T), min(dat$npgo.regr, na.rm=T)))
abline(h=mean(dat$npgo.regr, na.rm=T))
abline(v=1989, lty=2)
mtext("l", adj=0.05, line=-1.4, cex=1.1)

par(mar=c(1.5,1.5,1.5,1))


# finally, the NPGO-EOF2 plots
lim <- range(npgo.eof2.regr1, npgo.eof2.regr2)
z <- rep(NA, ncol(SST))
# z[!land] <- npgo.pattern1  # replace elements NOT corresponding to land with loadings!
# z[!land] <- npgo.eof.r1  # replace elements NOT corresponding to land with loadings!
z[!land] <- npgo.eof2.regr1  # replace elements NOT corresponding to land with loadings!

z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=new.col, zlim=c(lim[1], -lim[1]),
           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))


contour(x.t,y.t,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("m", adj=0.05, line=-1.4, cex=1.1)

z <- rep(NA, ncol(SST))
# z[!land] <- npgo.pattern2  # replace elements NOT corresponding to land with loadings!
# z[!land] <- npgo.eof.r2
z[!land] <- npgo.eof2.regr2
z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=new.col, zlim=c(lim[1], -lim[1]),
           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))


contour(x.t,y.t,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("n", adj=0.05, line=-1.4, cex=1.1)

par(mar=c(1.5,3,1.5,1))

plot(plot.dat$dec.yr, plot.dat$PDO.NPGO.cor, type="l", xlab="", ylab="PDO-NPGO correlation", xlim=xlim, col=cb[6])
abline(h=0)
abline(v=1989, lty=2)
mtext("o", adj=0.05, line=-1.4, cex=1.1)
dev.off()

