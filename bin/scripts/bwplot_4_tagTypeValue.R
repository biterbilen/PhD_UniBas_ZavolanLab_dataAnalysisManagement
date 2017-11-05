source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))
library(latticeExtra)

cargs <- commandArgs(trailingOnly=T)
clog <- F
ylabel <- "MIRZA Target Quality";
if (length(cargs)>1 & cargs[2]>0) {
	clog <- T
}
if (length(cargs)>2) {
	ylabel <- cargs[3]
}

a <- read.table(cargs[1],header=T,sep="\t")

plot.new(); pdf(paste(cargs[1],'.pdf',sep=''),h=8,w=8);
bwplot(value ~ tag | type, data=a, 
	scales=list(y=list(log=clog),x=list(rot=90)),
	yscale.components = yscale.components.log10ticks,
	strip = strip.custom(style=1,bg="white"),
	par.strip.text = list(cex=article.cex,font=article.font),
	ylab=ylabel,
	par.settings = article.theme, cex=article.cex,lwd=article.lwd,
	pch = "|",
	panel = function(x, y, ...) {
		panel.grid(h = -1, v = -1)
		panel.bwplot(x, y, ...)
		meds <- tapply(y, x, median)
		ylocs <- seq_along(meds)
		print(meds)
		print(ylocs)
		panel.segments(ylocs - 1/4, meds,
			ylocs + 1/4, meds,
			lwd = article.lwd, col = "red")
	}
)

dev.off()
