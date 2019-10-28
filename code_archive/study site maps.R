library(maps)
library(mapdata)


png("plots/study map.png", 6,6, units="in", res=300)
par(mar=c(0.2, 0.2, 0.2, 0.2))
plot(1:10, 1:10, xlim=c(150,245), ylim=c(30,62), type="n", yaxt="n", xaxt="n", ylab="", xlab="")


map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,80),add=T, lwd=1)
dev.off()