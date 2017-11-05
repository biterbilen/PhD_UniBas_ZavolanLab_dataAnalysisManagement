cargs <- commandArgs();

library(lattice)
#library(latticeExtra) #for axis.grid

#cargs=c(1,1,1,"~bilebi00/_EWSR1/Analysis/reverse_cum.rep")
if (length(cargs) == 4) {
	cargs[5] <- "T2C count";
	cargs[6] <- "Positions";
	cargs[7] <- "Reverse Cumulative";
} else if (length(cargs) == 5) {
	cargs[6] <- "Positions";
	cargs[7] <- "Reverse Cumulative";
} else if (length(cargs) == 6) {
	cargs[7] <- "Reverse Cumulative";
}

a <- read.table(cargs[4], header=T, check.names=F)

outfile <- paste(cargs[4],".pdf", sep="");

plot.new(); pdf(outfile, width=8, height=8);
tit <- paste(cargs[7], "distribution of", cargs[5]);
xyplot(log2(numberOfPositions) ~ log2(minT2C) |  proteinName, groups=sampleId, data=a, 
	xlab=paste("log2(", cargs[5], ")", sep=""),
	ylab=paste("log2(", cargs[6], ")", sep=""),
	#axis = axis.grid, 
	auto.key=list(title=tit), type=c('l','g'))
	#auto.key=list(title=tit, columns=5), type=c('l','g'))
dev.off()
 

