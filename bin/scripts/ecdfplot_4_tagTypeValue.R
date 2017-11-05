source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))
library(latticeExtra)

cargs <- commandArgs(trailingOnly=T)
clog <- F
xlabel <- "MIRZA Target Quality";
if (length(cargs)>1 & cargs[2]>0) {
	clog <- T
}
if (length(cargs)>2) {
	xlabel <- cargs[3]
}
ylabel <- "Empirical cumulative distribution"
a <- read.table(cargs[1],header=T)

nms <- unique(a$tag)
for (i in nms) {
	n <- which(a$tag==i);
	l <- length(n)
	id <- paste(a$tag[n], "\t(n=",l,")",sep="");
	a$tag[n] <- id
}

plot.new(); pdf(paste(cargs[1],"_xyplot.pdf",sep=''),h=8,w=8);
ecdfplot(~value|type,data=a, 
	xlab=xlabel,ylab=ylabel,
	scales=list(x=list(log=clog)),
	xscale.components = xscale.components.log10ticks,
	strip = strip.custom(style=1,bg="white"),
	par.strip.text = list(cex=article.cex,font=article.font),
	par.settings = article.theme, cex=article.cex,lwd=article.lwd,
	auto.key=list(space="top",columns=1, cex=article.cex,font=article.font, fontfamily=article.fontfamily,lwd=article.lwd, col=article.colors,lines=F),
	panel = function(x, y, ...) {
		panel.grid(h = -1, v = -1)
		panel.ecdfplot(x, ...)
	},
	groups=tag)
dev.off()

