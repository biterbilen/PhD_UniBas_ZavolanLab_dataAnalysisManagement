cargs <- commandArgs(trailingOnly = T);
#cargs <- c("a","PosteriorSum")

a <- read.table(cargs[1], header=F, sep="\t");
names(a) <- c("gid", "Name", "Score","Count");
head(a)

library(latticeExtra)


plot.new(); pdf(paste(cargs[1], '.PosteriorAverage.pdf', sep=''), h=8, w=8);
histogram(~(Score/Count), data=a, n=50, main=cargs[1], ylab=paste(" Density -", levels(a$Name[1]), sep=''), xlab=paste("Average ", cargs[2], sep=''));
#densityplot(~(PosteriorSum/Count), data=a, plot.points=F, n=20, main=cargs[1], ylab=paste(" Density -", levels(a$Name[1]), sep=''));
#xyplot(PosteriorSum/Count ~ Distance | Name, data=a, type=c("g"), main = a$Name[1], auto.key=list(columns=1, corner=c(0,0.5)),
#	layout=c(2,1)) +
#	layer(panel.loess(x, y, span = 0.05, col="black"))
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
