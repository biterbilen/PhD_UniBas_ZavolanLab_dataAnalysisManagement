maplot <- function(file, rheader) {

library(affy)
d <- read.delim(file, stringsAsFactors=TRUE, row.names=rheader, header=TRUE);
dp=d+10e-5;

plot.new(); png(paste(file,'_pre-norm.png',sep=''), width=1024, height=1024);
mva.pairs(as.matrix(dp), labels=names(d), cex=2, main="Pre-norm Expression of RNA-seq Expressed Genes")
dev.off();

#ma.plot( rowMeans(log2(y)), log2(y[, 1])-log2(y[, 2]), cex=1 ) 
#title("Pre norm: CLIP Expression in 100mer Windows (363 v 364)")
#dev.copy2pdf(file="maplot_crudeNorm_363_364.pdf");

library(preprocessCore)

#do a quantile normalization
x <- normalize.quantiles(as.matrix(dp));

plot.new(); png(paste(file,'_post-norm.png',sep=''), width=1024, height=1024);
mva.pairs(as.matrix(x), labels=names(d), cex=2, main="Post-norm Expression of RNA-seq Expressed Genes")
dev.off();

write.table(x, sep="\t", col.names=colnames(d), row.names=rownames(d), quote=F,
	file=paste(file,'.qn', sep=''));

#ma.plot( rowMeans(log2(x)), log2(x[, 1])-log2(x[, 2]), cex=1 ) 
#title("Post norm: RNA-seq Expression in 100mer Windows (59 v 60)")
#dev.copy2pdf(file="maplot_postNorm_59_60.pdf");

}

maplot(commandArgs()[4]);
quit();
