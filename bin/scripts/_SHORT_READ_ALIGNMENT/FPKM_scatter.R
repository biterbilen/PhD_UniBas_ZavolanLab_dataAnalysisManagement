cargs <- commandArgs();

#cargs <- c(1,1,1,"../Scheiffele/Cuffdiff_trusted_RNAseq/genes.fpkm_tracking");
#cargs <- c(1,1,1,"Zavolan/isoforms.selected_RepAB");
#cargs <- c(1,1,1,"Zavolan/genes.selected_RepAB");
#cargs <- c(1,1,1,"Zavolan/isoforms.selected_just_estimation_RepAB");
#cargs <- c(1,1,1,"Zavolan/genes.selected_just_estimation_RepAB");
library(latticeExtra);

#TODO what about pseudocount
aa <- read.table(cargs[4],header=T);
sel=apply(as.matrix(cbind(aa$FPKM.1, aa$FPKM.2)), 1, function (x) { x[1] > 0 && x[2] > 0 });
a <- aa[sel,]

plot.new(); pdf(paste(cargs[4], '.FPKM.pdf', sep=''),w=8,h=8)
xyplot(log2(FPKM.1) ~ log2(FPKM.2) | paste(g1,g2, sep=" VS "), data=a, type=c("g","p"), main=cargs[4] ) + 
	layer(panel.abline(a=0,b=1, rot=T, pos=3, lwd=2,col.line=2)) +
	layer(panel.ablineq(lm(y ~ x), rot=T, pos=3, lwd=2)) +
	#layer(panel.loess(x,y, r.sq=T, r.sq=T, rot=T, pos=3, col.line=4, lwd=2)) +
	layer(panel.key(c(sprintf('R=%.2f N=%d', cor(x,y), length(x))),lwd=2,cex=2))
dev.off()

