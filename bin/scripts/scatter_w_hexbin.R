library(latticeExtra)
library(hexbin)
source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))

cargs <- c("INTRON.relxlinkPos.relRegionPos.xlinkScore","x","y")
cargs <- commandArgs(trailingOnly=T)

a <- read.table(cargs[1],header=T,sep="\t")
xlabel <- cargs[2]
ylabel <- cargs[3]
xindex <- as.numeric(cargs[4])
yindex <- as.numeric(cargs[5])
tag <- cargs[6]
names(a)[1] <- "name"
names(a)[xindex] <- "x"
names(a)[yindex] <- "y"

plot.new(); pdf(paste(cargs[1],".",tag,".pdf",sep=''),h=8,w=8)
xyplot(y ~ x, data=a, 
	xlab=xlabel, ylab=ylabel,
	par.settings = article.theme, cex=article.cex*1,
	xscale.components = xscale.components.log10ticks,
	yscale.components = yscale.components.log10ticks,
	scale=list(log=F),main=cargs[1], 
	panel = function(x, y, ...){
		panel.superpose 
		panel.grid(h = -1, v = -1)
		panel.hexbinplot(x,y, ...)
		panel.loess(x,y,...,col=article.colors[3],lwd=article.lwd)
		panel.abline(h=0.5, ..., col=article.colors[3],lty="dashed", lwd=article.lwd)
#		panel.xyplot(x,y,...,col=article.colors[1])
		panel.text(max(x),max(y),sprintf("N=%.0f",length(x)),cex=article.cex,fontfamily=article.fontfamily,font=article.font)
	}
	)
dev.off()

