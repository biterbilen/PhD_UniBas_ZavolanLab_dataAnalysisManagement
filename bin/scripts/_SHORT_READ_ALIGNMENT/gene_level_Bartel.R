fargs <- commandArgs();
library(latticeExtra)
#cargs <- c(1,1,1,"../Guo/Trx_level/mrna_mir1_12hrmrna_mock_12hr-mrna_mir1_12hrmrna_mock_12hr.seeds", "mrna_mir1_12hrVSmrna_mock_12hr");

a <- read.table(cargs[4], header=T)

tag=paste(cargs[5],'cufflinksFPKM.pdf',sep="_");
plot.new(); pdf(tag,w=8,h=8)
xyplot(log(value_1+1) ~ log(value_2+1), data = a, main=tag) +
layer(panel.abline(a=0,b=1, col.line=2, lwd=1)) +
layer(panel.key(c(sprintf('R=%.2f N=%d', cor(x,y), length(x))),lwd=2,cex=2))
dev.off();

tag=paste(cargs[5],'bartelRPKM.pdf', sep="_");
plot.new(); pdf(tag,w=8,h=8)
xyplot(log(bartelRPKM1+1) ~ log(bartelRPKM2+1), data = a, main=tag) +
layer(panel.abline(a=0,b=1, col.line=2, lwd=1)) +
layer(panel.key(c(sprintf('R=%.2f N=%d', cor(x,y), length(x))),lwd=2,cex=2))
dev.off();

tag=paste(cargs[5],'bartelRPKM_logRatio.pdf', sep="_");
plot.new(); pdf(tag,w=8,h=8)
bwplot((log(bartelRPKM1+1) - log(bartelRPKM2+1)) ~ seedType,
	data=a, scales=list(x=list(rot=45)), main=tag) +
layer(panel.refline(h=0, col.line=2, lwd=1)) +
layer(panel.key(c(sprintf('N=%d', length(x))),lwd=2,cex=2))
dev.off()

tag=paste(cargs[5],"cufflinksFPKM_logRatio.pdf", sep="_");
plot.new(); pdf(tag,w=8,h=8)
bwplot((log(value_1+1) - log(value_2+1)) ~ seedType,
	data=a, scales=list(x=list(rot=45)), main=tag) +
layer(panel.refline(h=0, col.line=2, lwd=1)) +
layer(panel.key(c(sprintf('N=%d', length(x))),lwd=2,cex=2))
dev.off()

tag=paste(cargs[5],'cufflinksFPKM_logRatio_expCut.pdf', sep="_");
plot.new(); pdf(tag,w=8,h=8)
bwplot((log(value_1+1) - log(value_2+1)) ~ seedType,
	subset=value_1 > 1 & value_2 > 1,
	data=a, scales=list(x=list(rot=45)), main=tag) +
layer(panel.refline(h=0, col.line=2, lwd=1)) +
layer(panel.key(c(sprintf('N=%d', length(x))),lwd=2,cex=2))
dev.off()

quit()

tag=paste(cargs[5],'cufflinksFPKM_testStat_expCut.pdf', sep="_");
plot.new(); pdf(tag,w=8,h=8)
bwplot(test_stat ~ seedType,
	subset=value_1 > 0 & value_2 > 0,
	data=a, scales=list(x=list(rot=45)), main=tag) +
layer(panel.refline(h=0, col.line=2, lwd=1)) +
layer(panel.key(c(sprintf('N=%d', length(x))),lwd=2,cex=2))
dev.off()

