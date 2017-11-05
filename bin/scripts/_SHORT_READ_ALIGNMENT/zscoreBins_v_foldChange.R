cargs <- commandArgs();
library(latticeExtra)

#cargs <- c(1,1,1,"../Shivendra/CLIP_comparison_trusted_superClusters_gene/mRNAseq_siCTRL-A_HeLaVSmRNAseq_siEWSR1_HeLa.aggregate");

a <- read.table(cargs[4], header=T);

#edcf

#plot.new(); pdf(paste(cargs[4],'.edcflike.pdf',sep=''), h=8, w=8);
##xyplot(cumsum(abs(medianlog2FC)) ~ rank | CLIPtag, data=a, subset=CLIPtag=="PC_DIS3L2_1vPC_CFIm25_1andPC_DIS3L2_1vPC_CFIm68_2")+
##xyplot(cumsum(abs(medianlog2FC)) ~ rank | CLIPtag, data=a, subset=CLIPtag=="687-688_ALL_T2C.50.sp_cls") +
#xyplot(cumsum(abs(medianlog2FC)) ~ rank | CLIPtag, data=a, subset=CLIPtag=="687-688_mRNA_T2C.50.sp_cls") +
#layer(panel.loess(x, y, ..., col = 'black')) +
#layer(panel.key(c(sprintf('N=%d R=%.2f', length(x),cor(x,y,method="pearson"))),lwd=2,cex=1))
#dev.off()

plot.new(); pdf(paste(cargs[4],'.FC.pdf',sep=''), h=8, w=8);
#bwplot(medianlog2FC ~ equal.count(zscore,5,overlap=0) | CLIPtag, data=a, 
#bwplot(abs(medianlog2FC) ~ equal.count(zscore,5,overlap=0) | CLIPtag, data=a, 
#bwplot(medianlog2FC ~ cut(a$rank,breaks=c(1,100,250,500,1000,max(a$rank)),overlap=0) | CLIPtag, data=a, xlab="CLIP rank",
bwplot(~medianlog2FC | CLIPstat * CLIPtag, data=a,
#	panel=panel.violin, 
	scales=list(x=list(rot=45)), main=cargs[4],
	par.strip.text = list(cex=0.6)) +
layer(panel.key(c(sprintf('N=%d', length(x))),lwd=2,cex=1)) +
layer(panel.refline(v=0,lwd=2,cex=1))
dev.off();


plot.new(); pdf(paste(cargs[4],'.testStat.pdf',sep=''), h=8, w=8);
#bwplot(abs(testStat) ~ equal.count(zscore,5,overlap=0) | CLIPtag, data=a, 
#bwplot(testStat ~ cut(a$rank,breaks=c(1,100,250,500,1000,max(a$rank)),overlap=0) | CLIPtag, data=a, xlab="CLIP rank",  
bwplot(~testStat | CLIPstat * CLIPtag, data=a, 
#	panel=panel.violin, 
	scales=list(x=list(rot=45)), main=cargs[4], subscripts = T,
	par.strip.text = list(cex=.6)) +
#TODO this line to process group layer(panel.key(c(sprintf('%s; N=%d; CLIP N=%d', trellis.par.get("strip.shingle")$col[which.given], length(x), length(zscore[subscripts]>0))),lwd=2,cex=1),data=a) +
layer(panel.key(c(sprintf('N=%d', length(x))),lwd=2,cex=1)) +
layer(panel.refline(v=0,lwd=2,cex=1))
dev.off();

#xyplot(abs(medianlog2FC) ~ zscore | cut(a$zscore,breaks=c(-20,-2,0,2,20),overlap=0), data=a, 
#	subset=CLIPtag=="PC_DIS3L2_1vPC_CFIm25_1andPC_DIS3L2_1vPC_CFIm68_2", main=cargs[4]) +
#xyplot(abs(medianlog2FC) ~ zscore | CLIPtag, data=a, 
#	main=cargs[4]) +
#layer(panel.key(c(sprintf('N=%d; Sp R=%.2f', length(x),cor(x,y,method="spearman"))),lwd=2,cex=1)) +
#layer(panel.refline(h=0, col.line="darkgray", lwd=1)) +
#layer(panel.refline(v=0, col.line="darkgray", lwd=1))


#bwplot(~medianlog2FC | CLIPtag, data=a, subset=zscore>2,
#	panel=panel.violin, scales=list(x=list(rot=45)), main=cargs[4],
#	par.strip.text = list(cex=.4)) +
#layer(panel.key(c(sprintf('N=%d', length(x))),lwd=2,cex=1))

