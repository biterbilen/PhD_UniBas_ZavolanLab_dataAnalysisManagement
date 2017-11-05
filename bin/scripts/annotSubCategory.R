library(latticeExtra)
library(reshape)
source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))


cargs <- commandArgs(trailingOnly=T)
#cargs <- c("/import/bc2/home/zavolan/bilebi00/_CLIP/Analysis/Project_DIS3L2/xlinkEnrichment/Selected/categories")
#cargs <- c("categories.all")

a <- read.table(cargs[1], header=T,row.names="tag",check.names=F)
subsetCut <- 30;
if (length(cargs)>1) {
	subsetCut <- as.numeric(cargs[2]);
}
ylab <- "Gene Category Enrichment"
if (length(cargs)>2) {
	ylab <- cargs[3];
}
#a$preMRNA <- NULL
b <- melt(t(a))
d <- cast(b, X1 ~ X2)
#get ids and 
ids <- gsub(".all.annot","",grep("all",names(d),value=T))
#assign the shuffled one the ratio of the original over shuffled count
head(d)
dim(d)
sms <- apply(d,2,sum)
for (id in ids) {
	tag <- paste(id,".annot",sep='')
	tag2 <- paste(id,".all.annot",sep='')

	for (i in 1:nrow(d)) {
#		d[i,tag2] <- pbeta(0.5, d[i,tag], d[i,tag2],log=T)
#		d[i,tag2] <- d[i,tag] / d[i,tag2]
#		d[i,tag2] <- (d[i,tag] / sms[tag]) / (d[i,tag2] / sms[tag2])
		m <- matrix(c(d[i,tag],sms[tag],d[i,tag2],sms[tag2]),nr=2);
		d[i,tag2] <- -log(fisher.test(m, alternative="g")$p.value);
	}
#	d[tag2] <- d[tag] / d[tag2]
}
head(d)
k <- NULL
k$X1 <- d$X1
k["sum"] <- (d["288.annot"] + d["288.annot"]) / 2
#k["sum"] <- (d["445.all.annot"] + d["446.all.annot"]) / 2
mdn <- median(k$sum)
k <- data.frame(k)

e <- melt(d)
f <- subset(e,grepl("all", e$X2))
f <- merge(f, k, by="X1")
f$X2 <- factor(f$X2)
f$X2 <- gsub(".all.annot","",f$X2)
f$X1 <- gsub("UTR","'UTR",f$X1)
rownames(f) <- NULL
head(f)
plot.new(); pdf(paste(cargs[1],'.pdf',sep=''),h=8,w=8)
#barchart((-value) ~ reorder(factor(X1), sum), data=f,
#	group=X2, ylab="Enrichment (-log(p-value))",
barchart((value) ~ reorder(factor(X1), -sum), data=f, subset=sum>subsetCut, drop.unused.levels=T,
#	group=X2, ylab="Positions in top 1000 / All Positions",
	group=X2, ylab=ylab,
	par.settings = article.theme, lwd=article.lwd, cex=article.cex, origin=0,
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=article.cex,font=article.font),
	auto.key = list(space="top",points=F, rectangles=F, col=article.colors,columns=length(unique(f$X2)), cex=article.cex,font=article.font, fontfamily=article.fontfamily),
	scales=list(y=list(log=T)),
	yscale.components = yscale.components.log10ticks,
	panel = function(x,y,...) {
		panel.grid(h = -1, v = -1)
		panel.barchart(x, y, col=article.colors,...)
#		panel.abline(h=mdn,lwd=article.lwd,label="reference",adj=c(0,-0.1),cex=1,fontfamily=article.fontfamily,at=0.1,col.text="gray",col="gray",lty="dashed")
	}
	)
dev.off()

quit()

d <- merge(a,b,by="gid")
plot.new(); pdf(paste(cargs[1], cargs[2],'.pdf',sep=''),h=8,w=8)
xyplot(preMRNA.y ~ preMRNA.x, data=d, 
	main=paste(cargs),
	ylab=cargs[2], xlab=cargs[1],main="preMRNA",
	scales=list(y=list(log=T), x=list(log=T)))

xyplot(INTRON.y ~ INTRON.x, data=d, 
	ylab=cargs[2], xlab=cargs[1],main="INTRON",
	scales=list(y=list(log=T), x=list(log=T)))

dev.off()

