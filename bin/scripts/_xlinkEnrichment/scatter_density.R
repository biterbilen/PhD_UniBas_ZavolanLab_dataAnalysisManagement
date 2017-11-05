cargs <- commandArgs(trailingOnly=T);
#cargs <- c("b.gtf","count");


library(latticeExtra)
library(hexbin)
library(reshape)
source(file.path("/import/bc2/home/zavolan/bilebi00/_ARE/scripts/","theme.R"))
source(file.path("/import/bc2/home/zavolan/bilebi00/_xlinkEnrichment/scripts/","read_gtfPlus.R"))

#cargs <- c(1,1,1,"hg19.clusters.transcripts.genes", "11","14","/import/bc2/home/zavolan/bilebi00/_ARE/Analysis/PAchange_wIG/IG_density");
a <- getDataFromGtfPlus(cargs[1])
b <- a[,grep(cargs[2],names(a))];
head(b)
d <- log2(subset(b,mut_count>0 & mut_count.1>0))
head(d)

plot.new(); pdf(paste(cargs[1],'.splom.pdf',sep=''), h=8, w=8);
splom(d, panel=panel.hexbinplot,
	varname.cex=1.5, axis.text.cex=1.5, 
	type=c("g","p"),
 	main=cargs[1],
	auto.key=list(corner=c(1,1), cex=1.5,font=2),
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=1.5,font=2),
	par.settings = poster.theme,
	lwd=2,cex=1.5, font=2,
	lower.panel = function(x, y, ...){
		panel.hexbinplot(x, y, ...)#,
		panel.abline(a=0,b=1, lwd=2, col = 'darkgray')#,
		panel.loess(x, y, lwd=2, col = 'black')#,
		panel.key(c(sprintf('N=%d\nR=%.2f', length(x),cor(x,y) )),col='black',points=F,cex=1.5,font=2)
	}, 
	)
dev.off()

#------
ma <- reshape(d, direction="long", 
	varying=list(names(d)), v.names="value",
	timevar="type",times=names(d))

plot.new(); pdf(paste(cargs[1],'.density.pdf', sep=''), h=8, w=8)
densityplot(~value, data=ma, group=type, plot.points=F, 
	xlab=cargs[2],
	varname.cex=1.5, axis.text.cex=1.5, 
	type=c("g","p"),
 	main=cargs[1],
	auto.key=list(corner=c(1,1), cex=1.5,font=2),
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=1.5,font=2),
	par.settings = poster.theme,
	lwd=2,cex=1.5, font=2,
	span=0.01)
dev.off();

