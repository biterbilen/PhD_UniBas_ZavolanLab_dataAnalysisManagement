library(latticeExtra)
source(file.path("/import/bc2/home/zavolan/bilebi00/_ARE/scripts/","theme.R"))

cargs <- commandArgs(trailingOnly=T);
#cargs <- c("xlink.counts","xlinkedGene.counts","gene.counts","DIS3L2");

xlink <- read.table(cargs[1], header=T)
xlinkedGene  <- read.table(cargs[2], header=T)
gene <- read.table(cargs[3], header=T)

b <- merge(gene,merge(xlink,xlinkedGene,by=c("type")))

#filters
if (length(cargs) > 5) {
#	b <- b[grep("pseudogene|processed_transcript|protein_coding",a$type,invert=T,perl=T),]
	b <- b[grep(cargs[6],a$type,perl=T),]
}
b <- b[b$geneCount > 50,]
b$type <- gsub("_","",b$type)
b$lab <- "";
b

lis <- ((b$xlinkCount/b$totalXlinkCount > 0.02) | (b$xlinkedGeneCount/b$geneCount > 0.02))
lis
b$lab[lis] <- b$type[lis];
b
b$xlinkCount/b$totalXlinkCount

plot.new(); pdf(paste(cargs[5], 'annotation.pdf',sep=''),h=8,w=8);
xyplot(100 * (xlinkCount/totalXlinkCount) ~ 100 * (xlinkedGeneCount/geneCount), 
	main=cargs[4],
	ylab="Crosslinked Positions (%)",
	xlab="Crosslinked Genes (%)",
	data=b,
	type=c("g","p"), 
	auto.key=list(corner=c(1,1), cex=article.cex,font=article.font),
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=article.cex,font=article.font),
	par.settings = article.theme,
	lwd=article.lwd,cex=article.cex,font=article.font,pch="x" ) +
layer(panel.text(..., labels=b$lab, cex=article.cex, font=article.font))

dev.off();

quit();

plot.new(); pdf(paste(cargs[3],'relAnnot.pdf',sep=''),h=8,w=8);
xyplot(log10((countOfpernucs.DIS3L2/totalCountOfpernucs.DIS3L2) / (countOfpernucs.HuR/totalCountOfpernucs.HuR)) ~ 100 * (countOfGenes.DIS3L2-countOfGenes.HuR)/countOfGenes, data=a, type=c("g","p")) +
layer(panel.text(..., labels=a$type, cex=1.0))
dev.off();

