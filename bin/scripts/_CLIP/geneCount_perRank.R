library(latticeExtra)
library(hexbin)
source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))
cargs <- commandArgs(trailingOnly=T)
#cargs <- c("GenePerRank/308.309.genePerRank", "title")

a <- read.table(gzfile(cargs[1]), header=T, sep="\t")
dim(a) 
a$total_gene_count <- a$f1_gene_count + a$f2_gene_count - a$common_gene_count;

p1 <- xyplot((common_gene_count/total_gene_count) ~ rank | type, data=a, #[a$total_gene_count%%10==0,], 
	layout=c(1,length(unique(a$type))), main=paste("Recovered Genes versus Rank", cargs[2],sep=' - '),
	ylab="Ratio of Common Genes over Recovered Genes", xlab="Rank of Crosslinked Position",
	auto.key = list(space="top",points=F, col=article.colors,columns=length(unique(a$tag)), cex=article.cex,font=article.font, fontfamily=article.fontfamily),
	group=tag, origin=0,
	ylim=c(0,1),
	scale=list(y=list(log=F), x=list(log=T,cex=article.cex)), 
	xscale.components = xscale.components.log10ticks,
	panel = function(x,y,...,col=col) {
		panel.grid(h = -1, v = -1)
		panel.abline(h=0.5, col="gray",...)
		panel.xyplot(x, y, ...)
#		panel.text(1,0.8,"wo_filters",col=col)
#		panel.xyplot(x, y, ...,col=col)
	},
	par.settings = article.theme, lwd=article.lwd, cex=article.cex, pch="+",
	strip = strip.custom(style=1,bg="white"),
	par.strip.text = list(cex=article.cex,font=article.font))

if (length(cargs) > 2 & cargs[3] == "pdf") {
	plot.new(); pdf(paste(cargs[1],'.pdf',sep=''),h=8,w=8)
	print(p1)
	dev.off()
} else {
	plot.new(); png(paste(cargs[1],'.png',sep=''), h=720,w=720)
	print(p1)
	dev.off()
}
