cargs <- commandArgs();
library(latticeExtra);
library(hexbin);

cargs <- c(1,1,1,"../Stepanka/for_clip");

a <- read.table(cargs[4], header=T);

plot.new(); pdf(paste(cargs[4], '.exp.pdf', sep=''), h=8, w=8);
xyplot(log(value_1+1) ~ log(value_2+1) | annot_type, data=a, panel=panel.hexbinplot,main=as.character(a$name[1])) + 
	layer(panel.abline(a=0,b=1, lwd=2, col = 'black')) +
	layer(panel.loess(x, y, lwd=2, col = 'red')) +
	layer(panel.key(c(sprintf('N=%d; R=%.2f', length(x),cor(x,y) )),lwd=2,cex=2))
dev.off();

plot.new(); pdf(paste(cargs[4], '.lnFoldChange-zscore.pdf', sep=''), h=8, w=8);
xyplot(ln.fold_change. ~ score| annot_type, data=a, panel=panel.hexbinplot, subset=status=="OK"&value_1>1&value_2>1,main=as.character(a$name[1])) + 
	layer(panel.abline(a=0,b=1, lwd=2, col = 'black')) +
	layer(panel.loess(x, y, lwd=2, col = 'red')) +
	layer(panel.key(c(sprintf('N=%d; R=%.2f', length(x),cor(x,y) )),lwd=2,cex=2))
dev.off();

plot.new(); pdf(paste(cargs[4], '.testStat-zscore.pdf', sep=''), h=8, w=8);
xyplot(test_stat ~ score| annot_type, data=a, panel=panel.hexbinplot, subset=status=="OK"&value_1>1&value_2>1,main=as.character(a$name[1])) + 
	layer(panel.abline(a=0,b=1, lwd=2, col = 'black')) +
	layer(panel.loess(x, y, lwd=2, col = 'red')) +
	layer(panel.key(c(sprintf('N=%d; R=%.2f', length(x),cor(x,y) )),lwd=2,cex=2))
dev.off();

plot.new(); pdf(paste(cargs[4], '.zscore.pdf', sep=''), h=8, w=8);
xyplot(score ~ score2 | annot_type, data=a, panel=panel.hexbinplot,main=as.character(a$name[1])) + 
	layer(panel.abline(a=0,b=1, lwd=2, col = 'black')) +
	layer(panel.loess(x, y, lwd=2, col = 'red')) +
	layer(panel.key(c(sprintf('N=%d; R=%.2f', length(x),cor(x,y) )),lwd=2,cex=2))
dev.off();

plot.new(); pdf(paste(cargs[4], '.bw.pdf', sep=''), h=8, w=8);
bwplot(score ~ significant | annot_type, data=a, panel=panel.violin, subset=score>0&status=="OK",main=as.character(a$name[1])) +
	layer(panel.key(c(sprintf('N=%d', length(x) )),lwd=2,cex=2))
dev.off();


