cargs <- commandArgs()

#cargs <- c(1,1,1,"../../_ARE/Analysis/Reproducibility/combined.freq5mer", "5mer");

aa <- read.table(cargs[4], header=T, check.names=F);
a <- aa[,3:dim(aa)[2]];
a <- replace(a, a==0, NA);
cor.a <- cor(a, use = "pairwise.complete.obs");

write.table(cor.a,quote=F,file=paste(cargs[4],'.cor', sep=''),sep="\t",row.names=T, col.names=NA)

ord <- order.dendrogram(as.dendrogram(hclust(dist(cor.a))))

x <- cor.a
#x  <- t(as.matrix(scale(cor.a)))
dd.row <- as.dendrogram(hclust(dist(x)))
row.ord <- order.dendrogram(dd.row)

dd.col <- as.dendrogram(hclust(dist(t(x))))
col.ord <- order.dendrogram(dd.col)

library(latticeExtra)

ltheme <- canonical.theme(color = FALSE)      ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme, lwd=4)      ## set as default

outputfile <- paste(cargs[4],".pdf",sep="");
plot.new(); pdf(outputfile, h=8,w=8);
tit=paste("Correlation of ",cargs[5],"-freqs around non-reproduced T2C sites",sep="");
levelplot(x[row.ord, col.ord],
	par.strip.text = list(cex=1.5,font=2),
	par.settings=list(axis.text=list(font=2,cex=1.5),
		par.ylab.text=list(font=2,cex=1.5),
		par.xlab.text=list(font=2,cex=1.5),
		par.main.text=list(font=2,cex=1.5)),
	lwd=2,cex=1.5,
	aspect = "fill",
	xlab=list(label="samples"),
	ylab=list(label="samples"),
	main=list(label=tit),
#	col.regions=colorRampPalette(c("red","white","blue")),
	scales = list(x = list(rot = 90)),
	colorkey = list(space = "top"),
	legend =
	list(right =
		list(fun = dendrogramGrob,
			args =
			list(x = dd.col, ord = col.ord,
				side = "right",
				size = 10))))
dev.off() 

