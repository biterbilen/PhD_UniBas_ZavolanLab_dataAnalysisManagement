library(latticeExtra)
source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))

cargs <- c("1628_1637.5p_3p","1")
cargs <- commandArgs(trailingOnly=T)

a <- read.table(cargs[1],header=T,sep="\t");
minvalue <- as.numeric(cargs[2])
head(a)

plot.new(); pdf(paste(cargs[1],".pdf",sep=''),h=8,w=8)
xyplot(log2( X5p_count / X3p_count) ~ log2(X5p_count * X3p_count), data=a, subset=X5p_count>minvalue&X3p_count>minvalue,
	par.settings = article.theme, cex=article.cex*2,
	xscale.components = xscale.components.log10ticks,
	yscale.components = yscale.components.log10ticks,
	scale=list(log=F),main=cargs[1], 
	panel = function(x, y, ...){
		panel.superpose 
		panel.grid(h = -1, v = -1)
		panel.loess(x,y,...,col=article.colors[3],lwd=article.lwd)
		panel.xyplot(x,y,...,col=article.colors[1])
		panel.abline(h=0, ..., col=article.colors[3],lty="dashed", lwd=article.lwd)
	})
dev.off()


