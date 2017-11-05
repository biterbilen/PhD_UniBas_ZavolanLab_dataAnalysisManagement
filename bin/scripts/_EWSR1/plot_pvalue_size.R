library(lattice)
cargs <- commandArgs();

#cargs <- c(1,1,1,"a");

a <- read.table(cargs[4], header=T, sep="\t")
head(a,1)

plot.new(); pdf(paste(cargs[4],'.pdf',sep=''), h=8, w=8);
xyplot(-log10(pvalue) ~ log10(TFclusterCount), data=a,
	main="p-value v data size", type=c("g","p"))
dev.off();

