cargs <- c("288_289.distances.gz")
cargs <- commandArgs(trailingOnly=T)
source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))

#Ugly code: 2 repetitions
#1.
a <- read.table(gzfile(cargs[1]),header=F)
slop <- (a[1,"V3"] - a[1,"V2"] - 1 ) / 2;
a$distance <- a$V2+slop-a$V8
a[a$distance < 0,"distance"] <-  a[a$distance < 0,"distance"] * -1

b <- matrix(0*1:(slop*3+3),nrow=slop+1,ncol=3)
for (d in 1:(slop+1)) {
	x <- a[a$distance==d-1,"V5"]
	y <- a[a$distance==d-1,"V11"]
	if (length(x) > 30) {
		b[d,1] <- cor(x,y,method="pearson")
		b[d,2] <- length(x)
	} else {
		b[d,1] <- NA
		b[d,2] <- NA
	}
	b[d,3] <- d-1
}
a1 <- as.data.frame(b)
a1$type <- cargs[3]
#-----------------
#2.
a <- read.table(gzfile(cargs[2]),header=F)
slop <- (a[1,"V3"] - a[1,"V2"] - 1 ) / 2;
a$distance <- a$V2+slop-a$V8
a[a$distance < 0,"distance"] <-  a[a$distance < 0,"distance"] * -1

b <- matrix(0*1:(slop*3+3),nrow=slop+1,ncol=3)
for (d in 1:(slop+1)) {
	x <- a[a$distance==d-1,"V5"]
	y <- a[a$distance==d-1,"V11"]
	if (length(x) > 30) {
		b[d,1] <- cor(x,y,method="pearson")
		b[d,2] <- length(x)
	} else {
		b[d,1] <- NA
		b[d,2] <- NA
	}
	b[d,3] <- d-1
}
a2 <- as.data.frame(b)
a2$type <- cargs[4]

d <- rbind(a1,a2);
names(d) <- c("correlation", "N", "distance", "type")
tail(d)

library(latticeExtra)
plot.new(); pdf(paste(cargs[1],"combined.pdf",sep=""),h=8,w=8)
xyplot(correlation ~ distance, data=d,
	groups=type,
	xlab="Distance",
	ylab="Correlation of Crosslink Scores",
	par.settings = article.theme, col=article.colors,pch="+", lwd=article.lwd,
	xscale.components = xscale.components.log10ticks,
	scales=list(x=list(log=T)),
	auto.key=list(space="top",points=F,columns=2,col=article.colors),
	panel = panel.superpose,
	panel.groups = function(x,y,col, ...) {
		panel.grid(h = -1, v = -1)
		panel.xyplot(x, y, ...)
		panel.loess(x, y, col=col,...)
		panel.text(min(x),median(y,na.rm=T),pos=3,sprintf("N=%.0f",length(!is.na(y))),cex=article.cex,fontfamily=article.fontfamily,font=article.font,col=col)
	}
	)
dev.off()
