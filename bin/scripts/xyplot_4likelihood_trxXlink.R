library(latticeExtra)
source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))


cargs <- c("b")
cargs <- commandArgs(trailingOnly=T)

b <- read.table(cargs[1],header=F,sep="\t")
names(b) <- c("value","gid", "trx_xlinkScore")

plot.new(); pdf(paste(cargs[1],".pdf",sep=""), h=8,w=8)
xyplot(value ~ -trx_xlinkScore, data=b,
	main=cargs[2], ylab=cargs[3],xlab=cargs[4],
	par.settings = article.theme, cex=article.cex,
	xscale.components = xscale.components.log10ticks,
	yscale.components = yscale.components.log10ticks,
	pch=19,#scale=list(x=list(log=T)),
	panel = function(x,y,col,...) {
		panel.grid(h=-1,v=-1)
		panel.superpose
		panel.xyplot(x, y, ...)
#		panel.loess(x, y, ..., col = 'black')
		panel.abline(x, y, ..., col = 'gray')
		panel.text(median(x),0.5,sprintf("N=%.0f R=%.6f",length(x),cor(x,y)),fontfamily=article.fontfamily,pos=4,adj=c(2,2))
	}
	)
dev.off();

