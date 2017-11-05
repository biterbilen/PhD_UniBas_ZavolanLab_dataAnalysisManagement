library(latticeExtra)
source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))


cargs <- c("stats_varLen")
cargs <- c("len_lt_45")
cargs <- commandArgs(trailingOnly=T)

a <- read.table(cargs[1],header=T)

plot.new(); pdf(paste(cargs[1],".pdf",sep=""),h=8,w=8);
barchart(nuc1+nuc2+nuc3+nuc4+nuc5+nuc6+nuc7+nuc8+nuc9+nuc10+nuc11+nuc12+nuc13+nuc14+nuc15+nuc16+nuc17+nuc18+nuc19+nuc20 ~ id | nucleotide, 
	data=a,par.settings = article.theme,
	ylab="Reads with nucX over total number of unmapped reads",
	yscale.components = yscale.components.log10ticks,
	auto.key=list(space="top",columns=4,col=article.colors,rectangles=F,points=F),horizontal=F, scales=list(y=list(log=T),x=list(rot=90)))
  panel = function(x,y,...) {
		panel.grid(h = -1, v = -1)
		panel.barchart(x, y, ...)}
dev.off()

