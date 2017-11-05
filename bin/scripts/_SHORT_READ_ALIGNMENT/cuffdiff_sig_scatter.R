cargs <- commandArgs();
library(hexbin)
library(latticeExtra)
#cargs <- c(1,1,1,"../Stepanka/CLIP_comparison_trusted_superClusters_gene/totRNAseq-HEK-0_2VStotRNAseq-HEK-0_gene.diff");
#cargs <- c(1,1,1,"../_Stepanka/CLIP_comparison/totRNAseq-HEK-0VStotRNAseq-HEK-oeDis3L2.gene_exp.diff");
#cargs <- c(1,1,1,"../Stepanka/CLIP_comparison_trusted/totRNAseq-HEK-0VStotRNAseq-HEK-oeDis3L2.gene_exp.diff");

header <- c("test_id","gene_id","gene","locus","sample_1","sample_2","status","exp_1","exp_2","ln(fold_change)","test_stat","p_value","q_value","significant","tss");
a <- read.table(cargs[4], header=F, col.names=header);

##TODO hexbin plot is buggy
#plot.new(); pdf(paste(cargs[4], ".primary.trx.pdf", sep=""), h=8,w=8);
#xyplot(log2(V8) ~ log(V9), data=a, panel=panel.hexbinplot, main=cargs[4], type=c("g","p")) +  
#	layer(panel.key(c(sprintf('N=%d R=%.2f', length(x),cor(x,y))),lwd=2,cex=1)) +
#	layer(panel.refline(a=0,b=1, col.line=2)) +
#	layer(panel.loess(x, y, ..., col = 'black')) +
#	as.layer(xyplot(log2(V8) ~ log(V9), data=a, subset=V14=="yes"))   
#dev.off();

#FC and TS correlation
#plot.new(); pdf(paste(cargs[4], ".FCvTS.pdf", sep=""), h=8,w=8);
#xyplot((-V10) ~ V11, data=a, main=cargs[4], type=c("g","p"), subset=V11!="nan") +  
#	layer(panel.key(c(sprintf('N=%d Spearman R=%.2f', length(x),cor(x,y,use="na.or.complete",method="spearman"))),lwd=2,cex=1)) +
#	layer(panel.refline(a=0,b=1, col.line=2)) +
#	layer(panel.loess(x, y, ..., col = 'black')) 
#dev.off();

plot.new(); pdf(paste(cargs[4], ".levelplot.primary.trx.pdf", sep=""), h=8,w=8);
#levelplot(V11 ~ log2(V8)+ log2(V9) | equal.count(V11,6,overlap=0), data = a,main=cargs[4],
#levelplot(V11 ~ log2(V8)+ log2(V9) | cut(V11,breaks=c(min(a$V11,na.rm=T),-2,0,2,max(a$V11,na.rm=T)),overlap=0), data = a,main=cargs[4],
levelplot(test_stat ~ log2(exp_1)+ log2(exp_2), data = a,main=cargs[4],
#	par.settings = theEconomist.theme(with.bg = TRUE),
	lwd=2,cex=1.5,font=2,
	par.settings = custom.theme.2,
	panel = panel.levelplot.points, type = c("p", 
		"g"), aspect = "iso", prepanel = prepanel.default.xyplot) +
	layer(panel.key(c(sprintf('N=%d\nR=%.2f', length(x),cor(x,y))),lwd=2,cex=1)) 
#	layer(panel.refline(a=0,b=1, col.line=2)) +
#	layer(panel.loess(x, y, ..., col = 'black')) 
	#as.layer(xyplot(log2(V9) ~ log(V8), data=a, subset=V14=="yes"))   
dev.off();

#plot.new(); pdf(paste(cargs[4], ".densityplot.primary.trx.pdf", sep=""), h=8,w=8);
#histogram(~V11, nint=100, main=cargs[4], data=a) 
#dev.off();

