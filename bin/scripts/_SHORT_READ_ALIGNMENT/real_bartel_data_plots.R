cargs <- commandArgs();
#cargs <- c(1,1,1,"../Guo/hsa-miR-1.seeds.isoforms");
#cargs <- c(1,1,1,"Guo/hsa-miR-155.seeds.isoforms");
#cargs <- c(1,1,1,"deneme_");
library(latticeExtra);

aa <- read.table(cargs[4],header=T);
sel=apply(as.matrix(cbind(aa$FPKM.1, aa$FPKM.2)), 1, function (x) { x[1] > 0 && x[2] > 0 });
a <- aa[sel,]
	
plot.new(); pdf(paste(cargs[4], '.logRatio.pdf', sep=''),w=8,h=8)
bwplot(logRatio ~ seedType | paste(g1,g2, sep=" VS "), 
	data=a, layout=c(2,1), scales=list(x=list(rot=45)), main=cargs[4]) +
	layer(panel.refline(h=0, col.line=2, lwd=1)) +
	layer(panel.key(c(sprintf('N=%d', length(x))),lwd=2,cex=2))
#bwplot(logRatio ~ seedType | paste(g1,g2), data=a)
dev.off()

plot.new(); pdf(paste(cargs[4], '.testStat.pdf', sep=''),w=8,h=8)
bwplot(testStat ~ seedType | paste(g1,g2, sep=" VS "), 
	data=a, layout=c(2,1), scales=list(x=list(rot=45)), main=cargs[4]) +
	layer(panel.refline(h=0, col.line=2, lwd=1)) +
	layer(panel.key(c(sprintf('N=%d', length(x))),lwd=2,cex=2))
#bwplot(logRatio ~ seedType | paste(g1,g2), data=a)
dev.off()

plot.new(); pdf(paste(cargs[4], '.FPKM.pdf', sep=''),w=8,h=8)
xyplot(log2(FPKM.1) ~ log2(FPKM.2) | paste(g1,g2, sep=" VS ") * seedType , data=a, type=c("g","p"), auto.key=list(columns=3),title='seedType') + 
	layer(panel.ablineq(lm(y ~ x), rot=T, pos=3, lwd=2)) +
	#layer(panel.loess(x,y, r.sq=T, r.sq=T, rot=T, pos=3, col.line=4, lwd=2)) +
	layer(panel.key(c(sprintf('N=%d', length(x))),lwd=2,cex=2))
dev.off()

