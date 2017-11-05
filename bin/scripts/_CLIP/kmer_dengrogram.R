cargs <- commandArgs(trailingOnly=T);
#cargs <- c("filepat","colpat","out");
#cargs <- c("logpbeta","out","3");

count <-as.numeric(cargs[3]);
a <- NULL; #read.table(cargs[4], header=T)
nms <- NULL;
#XXX
for (i in 1:count) {
	f <- cargs[i+3];
	nms[i] <- cargs[i+3+count];
	if (i > 1) {
		b <- read.table(f, header=T)
		a <- merge(a,b,by="kmer");
	} else {
		a <- read.table(f, header=T)
	}
}
b <- a[,grepl(cargs[1],names(a))]
colnames(b) <- nms;
row.names(b) <- gsub("kmer.","",a$kmer)
#bs <- colSums(b)
#d <- matrix(rep(bs,dim(b)[1]),ncol=dim(b)[2],byrow=T)
#a <- as.matrix(b/d)
a <- b;

#a <- replace(a, a==0, NA);
#cor.a <- cor(t(a), use = "pairwise.complete.obs", method="spearman");

#write.table(cor.a,quote=F,file=paste(cargs[3],'.cor', sep=''),sep="\t",row.names=T, col.names=NA)

kmern <- 20 
if ( dim(a)[1] > kmern ) {
	x <- t(a[order(rowSums(a),decreasing=F)[1:kmern],])
} else {
	x <- t(a)
}

#ord <- order.dendrogram(as.dendrogram(hclust(dist(cor.a))))
ord <- order.dendrogram(as.dendrogram(hclust(dist(x))))
#x <- cor.a
#x  <- t(as.matrix(scale(cor.a)))
dd.row <- as.dendrogram(hclust(dist(x)))
row.ord <- order.dendrogram(dd.row)

dd.col <- as.dendrogram(hclust(dist(t(x))))
col.ord <- order.dendrogram(dd.col)

library(latticeExtra)

#ltheme <- canonical.theme(color = FALSE)      ## in-built B&W theme
#ltheme$strip.background$col <- "transparent" ## change strip bg
#lattice.options(default.theme = ltheme, lwd=4)      ## set as default

outputfile <- paste(cargs[2],".heatmap.pdf",sep="");
plot.new(); pdf(outputfile, h=8,w=8);
#tit=paste("Correlation of ",cargs[5],"-freqs around non-reproduced T2C sites",sep="");
tit<-cargs[2]
levelplot(x[row.ord, col.ord],
	par.strip.text = list(cex=1.5,font=2),
	par.settings=list(axis.text=list(font=2,cex=1.5),
		par.ylab.text=list(font=2,cex=1.5),
		par.xlab.text=list(font=2,cex=1.5),
		par.main.text=list(font=2,cex=1.5)),
	lwd=2,cex=1.5, font.family="mono",
	aspect = "fill",
	xlab=list(label=""),
	ylab=list(label=""),
	main=list(label=tit), 
	col.regions=colorRampPalette(c("red","white","blue")),
	scales = list(x = list(rot = 90)),
	colorkey = list(space = "left"),
	legend =
	list(right =
		list(fun = dendrogramGrob,
			args =
			list(x = dd.col, ord = col.ord,
				side = "right",
				size = 10)),
		top = 
		list(fun = dendrogramGrob,
			args =
			list(x = dd.row, ord = row.ord,
				side = "top", type="triangle",
				size = 2))))
dev.off()

