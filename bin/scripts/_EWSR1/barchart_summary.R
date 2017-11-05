cargs <- commandArgs(trailingOnly=T);
#cargs <- c("position.stats","1000", "0");
topN <- 100;
barchart.stack <- F;
if (length(cargs) > 1) {
	topN <- as.numeric(cargs[2])
}
if (length(cargs) > 2) {
	barchart.stack <- as.logical(as.numeric(cargs[3]))
}

a <- read.table(cargs[1], header=T, sep="\t");

if (length(unique(a$type)) > topN) {
	library(plyr)
	#sum library types
	d <- ddply(a,.(type),function(x) { data.frame(value=sum(x$value)) } );
	#get topN of the types
	types <- NULL
	types$type <- as.character(d[sort(d$value, decreasing=T,index.return=T)$ix,]$type[1:topN])
	#reduce a to types
	b <- merge(types,a);
	#calculate other percentage
	d <- ddply(b,.(lib),function(x) { data.frame(type="other",N=x$N,value=100-sum(x$value)) } );
	a <- rbind(b,d)
}

source(file.path("/import/bc2/home/zavolan/bilebi00/_ARE/scripts/","theme.R"))
#source(file.path("./","theme.R"))
library(latticeExtra)
plot.new(); pdf(paste(cargs[1], '.barchart.pdf', sep=''), h=8, w=8, colormodel="cmyk");
if (barchart.stack) {
	barchart(lib ~ value, groups=reorder(factor(type),-value), data=a, 
		stack = barchart.stack, auto.key=list(space="top", adj=1, cex=article.cex, font=article.font),
		strip = strip.custom(style=1,bg="gray"),
		par.strip.text = list(cex=article.cex, font=article.font),
		par.settings = article.theme,
		lwd=article.lwd,cex=article.cex, font=article.font,
		xlab="Proportion (%)")
} else {
	barchart(value ~ reorder(factor(type),-value), group=lib, data=a, stack = barchart.stack,
		par.settings = article.theme,
		origin=0,
		auto.key=list(corner=c(.9,.9),size=2,cex=article.cex, font=article.font),
		lwd=article.lwd, cex=article.cex, font=article.font,
		ylab="Proportion (%)", scales=list(x=list(rot=30),axs="i"))  + 
			layer(panel.grid(h=-1,v=0),under=T, theme=article.theme)
}
dev.off();

