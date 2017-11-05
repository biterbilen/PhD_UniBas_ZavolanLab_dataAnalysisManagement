cargs <- commandArgs(trailingOnly=T)
#cargs <- c("a.dist.hist");

a <- read.table(cargs[1],header=F);
names(a) <- c("count","x","distance")
head(a)
library(lattice)
plot.new(); pdf(paste(cargs[1],'.pdf',sep=''),h=4,w=4);
xyplot(100*count/max(a$count)~distance,data=a, scales=list(y=list(log=T)),type=c("g","h"),subset=distance>-100&distance<100,
	col="black"
	)
dev.off()

