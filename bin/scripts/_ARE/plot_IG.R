library(latticeExtra)
library(plotrix)
library(hexbin)
library(reshape)

cargs <- commandArgs();

print(length(cargs));
print(cargs[8:length(cargs)])

#cargs <- c(1,1,1,"hg19.clusters.transcripts.genes", "11","14","/import/bc2/home/zavolan/bilebi00/_ARE/Analysis/PAchange_wIG/IG_density");
a <- read.table(cargs[4], header=F);

ps <- 1;
indices <- as.numeric(cargs[5]):as.numeric(cargs[6]);

b <- a[,indices];
names(b) <- cargs[8:length(cargs)]

ltheme <- canonical.theme(color = FALSE)      ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme, lwd=4)      ## set as default

plot.new(); pdf(paste(cargs[7],'.splom.pdf',sep=''), h=8, w=8);
splom(log2(b+ps), panel=panel.hexbinplot,main=cargs[7], varname.cex=.8, axis.text.cex=.8, type=c("g","p"),
	par.strip.text = list(cex=1.5,font=2),
	pair.settings=list(axis.text=list(font=2,cex=1.5),
		par.ylab.text=list(font=2,cex=1.5),
		par.xlab.text=list(font=2,cex=1.5),
		par.main.text=list(font=2,cex=1.5)),
	lwd=2,cex=1,font=2,
	lower.panel = function(x, y, ...){
		panel.hexbinplot(x, y, ...)#,
		panel.abline(a=0,b=1, ..., col = 'darkgray')#,
		panel.loess(x, y, ..., col = 'black')#,
		panel.key(c(sprintf('N=%d\nR=%.2f', length(x),cor(x,y) )),...,col='black')
	}, 
	)
dev.off()

#------
a <- read.table(cargs[7], header=T);
m <- as.data.frame(a[,-1])
wmeans <- apply(m,2,mean)
wse <- apply(m,2,sd)/sqrt(dim(m)[1]);
## note that barplot() returns the midpoints of the bars, which plotCI
##  uses as x-coordinates
plot.new(); pdf(paste(cargs[7],'.CIplot.pdf',sep=''), h=8, w=8);
plotCI(barplot(wmeans,col="gray",ylim=c(0,max(wmeans+wse))),
		wmeans,wse,add=T, main=basename(cargs[7]))
title(main=basename(cargs[7]), ylab="IG");
dev.off();

ma <- reshape(a, direction="long", varying=list(names(a)[-1]), v.names="IG", 
	idvar=c("gid"), timevar="samplePair", times=names(a)[-1])

plot.new(); pdf(paste(cargs[7],'.density.pdf', sep=''), h=8, w=8)
densityplot(~IG | samplePair, data=ma, group=samplePair, plot.points=F, 
	par.strip.text = list(cex=1.5,font=2),
	par.settings=list(axis.text=list(font=2,cex=1.5),
		par.ylab.text=list(font=2,cex=1.5),
		par.xlab.text=list(font=2,cex=1.5),
		par.main.text=list(font=2,cex=1.5)),
	lwd=2,cex=1,font=2,
	na.rm=T, auto.key=list(space="top"), main=basename(cargs[7]))
dev.off();

