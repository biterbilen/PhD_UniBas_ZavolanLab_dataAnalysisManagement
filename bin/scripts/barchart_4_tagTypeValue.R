library(latticeExtra)
source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))

cargs <- commandArgs(trailingOnly=T)
#cargs <- c("288.closest_distances", "238_239")
#cargs <- c("288")

a <- read.table(cargs[1],header=T, stringsAsFactors=FALSE,sep="\t")
ylab <- "Ratio";
if (length(cargs) > 1) ylab <- cargs[2];
print(length(cargs))
print (ylab)

plot.new(); pdf(paste(cargs[1],'.pdf',sep=''),h=8,w=8)
barchart(value ~ reorder(factor(tag), -value)|type, data=a, drop.unused.levels=T,
	scales=list(x=list(rot=0)), stack=T,
	ylab=ylab,xlab="Group",group=tag,
	par.settings = article.theme, origin=0,
	strip = strip.custom(style=1,bg="white"),par.strip.text = list(cex=article.cex,font=article.font),
	auto.key = list(space="top",points=F, rectangles=F, col=article.colors,columns=1, cex=article.cex,font=article.font, fontfamily=article.fontfamily),
	prepanel = function(x, y,...) {
		list(xlim=levels(reorder(x, -y)))
	},
	panel = function(x,y,...) {
		panel.grid(h = -1, v = -1)
		panel.barchart(x, y, ...)
	}
	)

dev.off()


