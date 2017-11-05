library(latticeExtra)
source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))

cargs <- c("a")
cargs <- commandArgs(trailingOnly=T)
a <- read.table(cargs[1])

plot.new(); pdf(paste(cargs[1],".pdf",sep=""),h=8,w=8)
xyplot(V4 ~ V3|V1, data=a, groups=V2, 
	ylab="Reads per total", xlab="Position",
	par.settings = article.theme, cex=article.cex,
	auto.key = list (space="top",col=article.colors,points=F,drop.unused=T), col=article.colors, #alpha=article.alpha,
	scales = list(y=list(log=T)), yscale.components = yscale.components.log10ticks,
	strip = strip.custom(style=1,bg="white"),par.strip.text = list(cex=article.cex,font=article.font),
	type=c("g","s"),#npoints=3,
	panel= function(x,y,...) {
		panel.abline(v=0,lwd=article.lwd)	
		panel.abline(v=100,lty="dashed",lwd=article.lwd)	
		panel.xyplot(x,y,...)
	})
dev.off()

