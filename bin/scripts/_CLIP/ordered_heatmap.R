library(latticeExtra)
source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))

#given a table (w header and w rowheader), 
#clusters the rows and/or columns with a method (multidimentional euclidean distance/correlation)
cargs <- c("a")
cargs <- c("result")
cargs <- c("288_289.annot");
cargs <- commandArgs(trailingOnly=T)

#default arguments
file <- cargs[1]
op <- "none" #cor minkowski euclidean complete ...
clog <- 0
clusterCols <- 1
clusterRows <- 1
plotDendCols <- 1
plotDendRows <- 1
ps <- 1
drawxtick <- T
drawytick <- F
rn <- "gid"
op2 <- "complete" #clustering method: complete single average

#read arguments
if (length(cargs) > 1)  op          <- cargs[2]
if (length(cargs) > 2)  clog        <- as.numeric(cargs[3])
if (length(cargs) > 3)  clusterCols <- as.numeric(cargs[4])
if (length(cargs) > 4)  clusterRows <- as.numeric(cargs[5])
if (length(cargs) > 5)  plotDendCols<- as.numeric(cargs[6])
if (length(cargs) > 6)  plotDendRows<- as.numeric(cargs[7])
if (length(cargs) > 7)  ps          <- as.numeric(cargs[8]) 
if (length(cargs) > 8)  drawxtick   <- as.logical(cargs[9])
if (length(cargs) > 9)  drawytick   <- as.logical(cargs[10])
if (length(cargs) > 10) rn          <- cargs[11]
if (length(cargs) > 11) op2         <- cargs[12]

#read file
aa <- read.table(file,header=T,row.names=rn,check.names=F,sep="\t");
head(aa)
#aa <- aa[grep("[KB]",rownames(aa),perl=T),]
#print(dim(aa))

#convert log
if (clog > 0) {
	a <- log2(as.matrix(aa)+ps)
	tag = " (log2)"
} else {
	a <- as.matrix(aa)
	tag = ""
}
head(a,2)

p=length(names(aa))
#calculate pairwise distances and dendrogram
if ( op != "none" & clusterCols > 0) {
	if ( op == "cor" ) { 
		cor.c <- cor(a, use = "pairwise.complete.obs", method="pearson");
		dd.c <- as.dendrogram(hclust(dist(cor.c)))
	} else {
		dd.c <- as.dendrogram(hclust(dist(t(a), method=op,p=p),method=op2))
	}
	dd.c.ord <- order.dendrogram(dd.c)
} else {
	dd.c.ord <- 1:dim(a)[2]
}
if (op != "none" & clusterRows > 0) {
	if ( op == "cor" ) {
		cor.r <- cor(t(a), use = "pairwise.complete.obs", method="pearson");
		dd.r <- as.dendrogram(hclust(dist(cor.r)))
	} else {
		dd.r <- as.dendrogram(hclust(dist(a, method=op,p=p),method=op2))
	}
	dd.r.ord <- order.dendrogram(dd.r)
	#TODO WRITE???
	#		write.table(cor.r[dd.r.ord,],quote=F,file=paste(file,'.cor', sep=''),sep="\t",row.names=T)
} else {
	dd.r.ord <- 1:dim(a)[1]
}

#prep list of dendrogram
#columns
dend = list()
if (plotDendCols > 0) {
	dd.c.list =  list(fun = dendrogramGrob,
		args =
		list(x = dd.c, ord = dd.c.ord,
			cex=0.5,
			side = "top",
			size = 2))
	dend <- list(top=dd.c.list)
}
#rows
if (plotDendRows > 0) {
	dd.r.list = list(fun = dendrogramGrob,
		args =
		list(x = dd.r, ord = dd.r.ord,
			cex=0.5,
			side = "right",
			size = 15))
	dend <- list(right=dd.r.list)
} 
if (plotDendRows == 1 & plotDendCols == 1) {
	dend <- list(top=dd.c.list,right=dd.r.list)
}

h <- 8
w <- 8
if(length(dd.r.ord) > 100) { w <- 2 }

print (c(h,w))
p <- levelplot(t(a[dd.r.ord,dd.c.ord]),
	par.settings = article.theme, lwd=article.lwd, cex=article.cex,
	strip = strip.custom(style=1,bg="white"),
	par.strip.text = list(cex=article.cex,font=article.font),
	aspect = "fill",
	xlab=list(label=""),
	ylab=list(label=""),
	sub=list(label=paste(file,tag,"\n(N=",dim(a)[1],")",sep=" "),pos=0),
	#col.regions=colorRampPalette(c("white","darkred")),
	col.regions=colorRampPalette(c("white","black")),
	scales = list(x = list(rot = 90,draw=drawxtick), y=list(draw=drawytick)),
	colorkey = list(space = "bottom"),
	legend = dend
	)

outputfile <- paste(cargs[1],".dendrogram.pdf",sep="");
plot.new(); pdf(outputfile, h=h,w=w);
plot(p)
dev.off() 


