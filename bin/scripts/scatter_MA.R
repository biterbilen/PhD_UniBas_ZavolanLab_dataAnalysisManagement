library(latticeExtra)
library(hexbin)
source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))

cargs <- c("phastCons","0")
cargs <- commandArgs(trailingOnly=T)

a <- read.table(cargs[1],header=T,sep="\t");
minvalue <- as.numeric(cargs[2])
xlabel <- "Average";
ylabel <- "Ratio (log2)";
names(a) <- c("id","x","y");
head(a)

b <- subset(a,x>minvalue&y>minvalue)
plot.new(); pdf(paste(cargs[1],".pdf",sep=''),h=8,w=8)
hexbinplot(log2(b$x / b$ y) ~ 0.5 * log2(b$x * b$y), 
	xlab=xlabel,ylab=ylabel,
#	par.settings = article.theme, #cex=article.cex, pch=".",
	main=cargs[1], style="nested.centroids",border=F,
	panel = function(x, y, ...){
		panel.superpose 
		panel.grid(h = -1, v = -1)
		#panel.xyplot(x,y,...,col=article.colors[1])
		panel.hexbinplot(x,y,...)
		panel.loess(x,y,...,col=article.colors[1],lwd=article.lwd*2)
		panel.abline(h=median(y), ..., col=article.colors[1],lwd=article.lwd*2)
		panel.abline(v=median(x), ..., col=article.colors[1],lwd=article.lwd*2)
		panel.abline(h=0, ..., col=article.colors[1],lty="dashed", lwd=article.lwd*2)
		panel.text(max(x),max(y),sprintf("N=%.0f",length(x)),cex=article.cex,fontfamily=article.fontfamily,font=article.font)
	})
dev.off()


