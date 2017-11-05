cargs <- c(".distance")
cargs <- commandArgs(trailingOnly=T)
source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))

a <- read.table(gzfile(cargs[1]),header=F)

distances <- unique(a$distance) 
b <- NULL
d <- distances[1]
for (d in distances) {
	x <- a[a$distance == d,]
	if (size(x)[1] > 100) {
		b <- c(d,cor(x[2],x[3],method="pearson"))
	} else {
		b[d,1] <- NA
		b[d,2] <- NA
	}
	b[d,3] <- d-1
}
a1 <- as.data.frame(b)

d <- rbind(a1);
names(d) <- c("correlation", "N", "distance", "type")
tail(d)

library(latticeExtra)
plot.new(); pdf(paste(cargs[1],".autocorrelation.pdf",sep=""),h=8,w=8)
xyplot(correlation ~ distance, data=d,
	groups=type,
	xlab="Distance",
	ylab="Correlation of Crosslink Scores",
	par.settings = article.theme, col=article.colors,pch=".", lwd=article.lwd,
	xscale.components = xscale.components.log10ticks,
	scales=list(x=list(log=T)), types=c("g","s"),
	auto.key=list(space="top",points=F,columns=2,col=article.colors),
	panel = panel.superpose,
	panel.groups = function(x,y,col, ...) {
		panel.grid(h = -1, v = -1)
		panel.xyplot(x, y, ...)
		#panel.loess(x, y, col=col,...)
		panel.text(min(x),median(y,na.rm=T),pos=3,sprintf("N=%.0f",length(!is.na(y))),cex=article.cex,fontfamily=article.fontfamily,font=article.font,col=col)
	}
	)
dev.off()
