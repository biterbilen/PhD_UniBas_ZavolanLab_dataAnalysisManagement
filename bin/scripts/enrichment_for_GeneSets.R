library(latticeExtra)
source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))

cargs <- c("a")
cargs <- c("TRANSCRIPT.288_289.kegg")
cargs <- commandArgs(trailingOnly=T)

a <- read.table(cargs[1],sep="\t")
head(a,1)

plot.new(); pdf(paste(cargs[1],'.pdf',sep=''),h=8,w=8)
barchart(reorder(V1, -log10(V8)) ~ -log10(V8), data=a, 
	par.settings = article.theme, col=article.colors[3],
	ylab="",xlab="Enrichment\n(-log10(p-value))",
	subset=V8<0.05,origin=0,
	panel = function(x,y,...) {
		panel.grid(h = -1, v = -1)
		panel.barchart(x, y, ...)
		panel.abline(h=0,v=2,lwd=article.lwd,lty="dashed",...,col.line="red")
	})
dev.off()

