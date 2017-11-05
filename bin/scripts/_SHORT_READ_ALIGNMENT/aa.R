cargs <- commandArgs();
library(latticeExtra)

cargs <- c(1,1,1,"Hela.aggregate");

a <- read.table(cargs[4], header=T);

plot.new(); pdf(paste(cargs[4],'.pdf',sep=''), h=8, w=8);
bwplot(medianlog2FC ~ factor(cut(a$zscore,5)) | tag, data=a, 
	groups=cut(a$zscore, 5), panel=panel.violin, scales=list(x=list(rot=45)), main=cargs[4])
dev.off();

