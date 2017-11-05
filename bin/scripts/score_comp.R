source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))
library(latticeExtra)
cargs <- commandArgs(trailingOnly=T)
#cargs <- c("all.1.bedplus.gz.nuc.gz");
#cargs <- c("mirRevCompSeedCount.stats", "Cumulative Reverse Complement Seed Count Rate")
 
a <- read.table(gzfile(cargs[1]), header=T,sep="\t")
head(a)
p1 <- xyplot(100 * rate ~ rank | as.factor(tag), data=a, #head(a,1000),
	layout=c(1,length(unique(a$tag))),
	par.settings = article.theme, lwd=article.lwd, cex=1,col=article.colors,
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=article.cex,font=article.font),
	scales=list(x=list("log"=T),y=list("log"=F)),
	xlab="Rank of Crosslinked Position", main=cargs[2], ylab="",
	auto.key = list(space="top",points=F, rectangles=F, col=article.colors,columns=length(unique(a$type)), cex=article.cex,font=article.font, fontfamily=article.fontfamily),
	panel = function(x,y,...) {
		panel.grid(h = -1, v = -1)
		panel.xyplot(x, y, ...)
	},
	groups=type, pch='o')

plot.new(); png(paste(cargs[1],'.png',sep=''), h=720, w=720); 
print(p1)
dev.off();

if (length(cargs) > 2 & cargs[3] == "pdf") {
	plot.new(); pdf(paste(cargs[1],'.pdf',sep=''), h=8, w=8); 
	print(p1)
	dev.off();
}
