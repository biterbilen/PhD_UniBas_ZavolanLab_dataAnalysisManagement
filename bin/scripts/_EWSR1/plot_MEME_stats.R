cargs <- commandArgs();
library(latticeExtra)

#cargs <- c(1,1,1,"200.a.stats");
#cargs <- c(1,1,1,"~bilebi00/_EWSR1/Analysis/MEME/stats");
a <- read.table(cargs[4], header=F);

plot.new(); pdf(paste(cargs[4],".pdf", sep=""), h=8, w=8);
i=1
print(densityplot((as.vector(a[-1,i])), lwd=3, col="darkgray", 
		xlim=c(0.9*min(a[,i]),1.1*max(a[,i])),main=cargs[4],
		xlab='MEME Motif Site Count') +
	layer(panel.refline(v=a[1,i], col.line="black", lwd=3)),
	split = c(1, 1, 1, 2)) 

i=2
print(densityplot(log10(as.vector(a[-1,i])), lwd=3, col="darkgray", 
		xlim=c(0.9*min(log10(a[,i])),1.1*log10(max(a[,i]))),main=cargs[4],
		xlab='MEME Motif log10 E-value') +
	layer(panel.refline(v=log10(a[1,i]), col.line="black", lwd=3)),
	split = c(1, 2, 1, 2), newpage = FALSE) 

dev.off();

