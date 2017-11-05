cargs <- commandArgs(trailingOnly=T)
#cargs <- "AGO2.posterior.stats"

a <- read.table(cargs[1],header=T) 
i <- a$posterior == as.numeric(cargs[2])
a$lab <- ""
a$lab[i] <- a$posterior[i]

library(latticeExtra)
source("/import/bc2/home/zavolan/bilebi00/_ARE/scripts/theme.R")

plot.new();pdf(paste(cargs[1],'.pdf',sep=''),h=8,w=8)
xyplot(common_over_union_percent ~ posterior, data=a,
	xlab="Posterior Cutoff",ylab="Reproduced Positions in Replicates (%)",
	main=paste("Positional Reproducibility of Posterior Cutoff for ", as.character(a$proteinname[1]), sep=""),
	par.settings = poster.theme,
	cex=1.5,font=2,lwd=2,col="black",
	scales=list(x=list(log=F)),type=c("g","p","l")) +
#	layer(panel.text(..., labels=a$lab, cex=1.5, font=2, srt=0, pos=3, col="black"))
	as.layer(xyplot(common_over_union_percent ~ lab, data=a, cex=5, pch=16, col="black"))
dev.off()

