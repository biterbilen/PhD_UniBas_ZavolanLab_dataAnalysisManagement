library(latticeExtra)
source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))


cargs <- c("a")
cargs <- commandArgs(trailingOnly=T)

a <- read.table(cargs[1], header=F,sep="\t")
names(a) <- c("group","distance","clusterCount")
a$clusterFraction <- a$clusterCount / a[1,3]
head(a)

plot.new(); pdf(paste(cargs[1],".pdf",sep=''),h=8,w=8)
xyplot(clusterFraction ~ distance, data=a, 
	groups=group,
	type=c("l"),pch="+",
	par.settings = article.theme, lwd=article.lwd, cex=article.cex,
	strip = strip.custom(style=1,bg="white"),
	par.strip.text = list(cex=article.cex,font=article.font),
	auto.key = list(space="top",points=F, rectangles=F, col=article.colors,cex=article.cex,font=article.font, fontfamily=article.fontfamily),
	ylab=cargs[2], xlab=cargs[3], main=cargs[4],
	xscale.components = xscale.components.log10ticks,
	yscale.components = yscale.components.log10ticks,
	panel = function (x,y,...){
		panel.grid(h = -1, v = -1)
		panel.xyplot(x,y,...)
	},
	scales=list(x=list(log=T),y=list(log=F))


	)
dev.off()

