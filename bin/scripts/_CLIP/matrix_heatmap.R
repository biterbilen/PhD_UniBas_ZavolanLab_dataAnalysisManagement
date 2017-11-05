library(latticeExtra)

cargs <- c("result")
cargs <- commandArgs(trailingOnly=T)
aa <- read.table(cargs[1],header=T,row.names="gid",check.names=F,sep="\t");

ps <- 1
a <- log2(as.matrix(t(aa))+ps)
#a <- as.matrix(t(aa))
head(a,2)
##a <- replace(a, a==0, NA);
#cor.a <- cor(a, use = "pairwise.complete.obs", method="pearson");
##cor.a <- cor(a, use = "pairwise.complete.obs", method="spearman");
#head(cor.a)
##write.table(cor.a,quote=F,file=paste(cargs[4],'.cor', sep=''),sep="\t",row.names=T, col.names=NA)

#ord <- order.dendrogram(as.dendrogram(hclust(dist(cor.a))))

#x <- cor.a
#x  <- t(as.matrix(scale(cor.a)))
#dd.row <- as.dendrogram(hclust(dist(x)))
#row.ord <- order.dendrogram(dd.row)

#dd.col <- as.dendrogram(hclust(dist(t(x))))
#col.ord <- order.dendrogram(dd.col)

#head(x)

#tit=paste("Correlation of ",cargs[5],"-freqs around non-reproduced T2C sites",sep="");
tit<-"";
#p <- levelplot(a[,row.ord ],
p <- levelplot(a,
	par.strip.text = list(cex=1.5,font=2),
	par.settings=list(axis.text=list(font=2,cex=1.0),
		par.ylab.text=list(font=2,cex=1.5),
		par.xlab.text=list(font=2,cex=1.5),
		par.main.text=list(font=2,cex=1.5)),
	lwd=2,cex=1.5,
	aspect = "fill",
	xlab=list(label=""),
	ylab=list(label=""),
	main=list(label=tit), 
	col.regions=colorRampPalette(c("white","black")),
	scales = list(x = list(rot = 90)),
	colorkey = list(space = "top")
#	legend =
#	list(right =
#		list(fun = dendrogramGrob,
#			args =
#			list(x = dd.col, ord = col.ord,
#				side = "right",
#				size = 10)))
	)

outputfile <- paste(cargs[1],".heatmap.pdf",sep="");
plot.new(); pdf(outputfile, h=8,w=8);
plot(p)
dev.off() 


