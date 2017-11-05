cargs <- commandArgs();

#cargs <- c(1,1,1,"array-cuffdiff_splicing_per_gene.stat");
a <- read.table(cargs[4], header=T);

library(latticeExtra)

of = paste(cargs[4], '.pdf', sep="");
plot.new(); pdf(of, h=8, w=8);
xyplot(log(maxsqrtJS) ~ log(maxabsIrank), data=a, subset=maxsqrtJS>0 & maxabsIrank>0, main=cargs[4]) +
#xyplot(maxsqrtJS ~ (abs(maxabsIrank)), data=a) +
	layer(panel.key(c(sprintf('N=%d; Spearman R=%.2f', length(x),cor(x,y,method="spearman") )),lwd=2,cex=2))
#splom(abs(a[2:7]));
dev.off()

