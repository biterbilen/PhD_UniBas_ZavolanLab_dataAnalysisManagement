library(latticeExtra)
source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))


cargs <- c("hha1")
cargs <- commandArgs(trailingOnly=T)

ylabel <- "G-nucleotide Count of 41-nucleotide-long Regions"
a <- read.table(cargs[1],header=T)
binSize <- as.numeric(cargs[2]);
if (length(cargs)>2) ylabel <- cargs[3];


p <- bwplot(value ~ floor(rank/binSize), data=a,
	xlab = paste("Sorted Bin (size=", binSize, ")", sep=""),
	ylab = ylabel,
	scales=list(y=list(log=F),rot=0),
	yscale.components = yscale.components.log10ticks,
	strip = strip.custom(style=1,bg="white"),
	par.strip.text = list(cex=article.cex,font=article.font),
	par.settings = article.theme, cex=article.cex,lwd=article.lwd,
	pch = "|",
	horizontal=F,
	panel = function(x, y, ...) {
		panel.superpose
		panel.grid(h = -1, v = -1)
		panel.bwplot(x, y, ...)
		meds <- tapply(y, x, median)
		print(meds)
		cnts <- tapply(y, x, length)
		print(cnts)
		ylocs <- seq_along(meds)
		panel.segments(ylocs - 1/4, meds,
			ylocs + 1/4, meds,
			lwd = article.lwd, col = "red")
	}
)
plot.new(); pdf(paste(cargs[1],".pdf",sep=""),h=6,w=8);
print(p)
dev.off()

