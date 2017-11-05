cargs <- commandArgs(trailingOnly=T)
#cargs <- "AGO2.overlap.stats"

a <- read.table(cargs[1],header=T) 
head(a)
a$sample1name <- gsub(".postcut.sorted.bed","", a$sample1name);
a$sample2name <- gsub(".postcut.sorted.bed","", a$sample2name);

library(latticeExtra)
source("/import/bc2/home/zavolan/bilebi00/_ARE/scripts/theme.R")
head(a)

plot.new();pdf(paste(cargs[1],'.pdf', sep=''),h=8,w=8)
barchart(intersection_over_topN_percent ~ as.factor(topN) | paste(sample1name,sample2name,sep=" and "),
	xlab="Top N",ylab="Reproduced Positions in Replicates (%)",
	layout=c(1,2), 
	main=paste("Positional Reproducibility of Top N for ", as.character(a$proteinname[1]), sep=""),
#barchart(intersection_over_topN_percent ~ paste(sample1name,sample2name,sep="\n") | as.factor(paste("topN=",topN,sep="")),
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=1.5,font=2),
	par.settings = poster.theme,
	cex=1.5,font=2, lwd=2,
	group=distance, data=a, 
	scales=list(y=list(limits=c(0,100))), 
	auto.key=list(space="top",title="distance",columns=length(unique(a$distance)),cex=1.5,font=2,cex.title=1.5,font.title=2)) +
	layer(panel.grid(h=-1,v=0))	
dev.off()


