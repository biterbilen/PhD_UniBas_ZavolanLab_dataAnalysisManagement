library(latticeExtra)
source(file.path("/import/bc2/home/zavolan/bilebi00/_ARE/scripts/","theme.R"))

cargs <- c("genecount.enrichment")
cargs <- commandArgs(trailingOnly=T)

a <- read.table(cargs[1],header=T,sep="\t",row.names="type")
a <- as.data.frame(t(a))
a
a$ratio <- a[,1] / a["ALL",1]
a$enrichment <- a$ratio
all1 <- a["ALL",1]
all2 <- a["ALL",2]
d <- a[1:(nrow(a)-1),]
for (i in 1:nrow(d)) {
	m <- matrix(c(d[i,1], all1-d[i,1], d[i,2], all2-d[i,2]),nr=2);
	d[i,"enrichment"] <- -log10(fisher.test(m, alternative="g")$p.value);
}
d

plot.new(); pdf(paste(cargs[1], '.pdf',sep=''),h=8,w=8);
xyplot(enrichment ~ (100*ratio), data=d,
	ylab="Gene Enrichment (-log10(p-value))\n(EWSR1 / HuR)",
	xlab="Fraction of Genes in Region (%)",
	type=c("g","p"),
	auto.key=list(corner=c(1,1), cex=article.cex,font=article.font),
	strip = strip.custom(style=1,bg="white"),
	par.strip.text = list(cex=article.cex,font=article.font),  
	par.settings = article.theme, col=article.colors,
	lwd=article.lwd,cex=article.cex,font=article.font,pch=19 ) +
layer(panel.text(..., labels=rownames(d), cex=article.cex, font=article.font,pos=3))

dev.off();

