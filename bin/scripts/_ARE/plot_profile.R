cargs <- commandArgs(trailingOnly = T);
#cargs <- c("profile")

a <- read.table(cargs[1], header=F);
names(a) <- c("Distance", "Name", "PosteriorSum","Count");
head(a)

library(latticeExtra)


plot.new(); pdf(paste(cargs[1], '.PosteriorAverage.pdf', sep=''), h=8, w=8);
xyplot(PosteriorSum/Count ~ Distance | Name, data=a, type=c("g","p"), main = a$Name[1], auto.key=list(columns=1, corner=c(0,0.5)),
	layout=c(1,2)) +
	layer(panel.loess(x, y, span = 0.05, col="black"))
dev.off()

quit();






plot.new(); pdf(paste(cargs[1], '.PosteriorSum.pdf', sep=''), h=8, w=8);
xyplot(PosteriorSum ~ Distance | Name, data=a, type=c("g"), main = a$Name[1], auto.key=list(columns=1, corner=c(0,0.5)),
	layout=c(2,1)) +
	layer(panel.loess(x, y, span = 0.05, col="black"))
dev.off()

plot.new(); pdf(paste(cargs[1], '.Count.pdf', sep=''), h=8, w=8);
xyplot(Count ~ Distance | Name, data=a, type=c("g"), main = a$Name[1], auto.key=list(columns=1, corner=c(0,0.5)),
	layout=c(2,1)) +
	layer(panel.loess(x, y, span = 0.05, col="black"))
dev.off()
