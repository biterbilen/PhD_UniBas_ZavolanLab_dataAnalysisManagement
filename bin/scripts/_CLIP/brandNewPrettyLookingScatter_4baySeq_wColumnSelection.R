library(latticeExtra)
source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))

cargs <- c("a","pernucExp","0.5","gid", "EWSR1")
cargs <- c("genes_DE",".pernucExp","0.5","gid", "EWSR1")
cargs <- commandArgs(trailingOnly=T)

likelihood <- 0.5
lblsn <- "gid"
labelName <- NULL;
if (length(cargs)>2) likelihood <- as.numeric(cargs[3])
if (length(cargs)>3) lblsn <- cargs[4]
if (length(cargs)>4) labelName <- cargs[5]

a <- read.table(cargs[1],header=T);
head(a)

if (length(grep("group",names(a))) == 0) {
	a$group <- "-"
	a$group[a$log2FC<0&a$Likelihood>likelihood] <- "upregulated"
	a$group[a$log2FC>0&a$Likelihood>likelihood] <- "downregulated"

	write.table(a[a$group=="upregulated",c("gid","ids")], file=paste("upregulated"),sep="\t",quote=F,row.names=F,col.names=F)
	write.table(a[a$group=="downregulated",c("gid","ids")], file=paste("downregulated"),sep="\t",quote=F,row.names=F,col.names=F)

	n <- length(a[a$group == "downregulated","group"])
	a$group[a$group=="downregulated"] <- paste("downregulated", " (n=",n,")",sep="")
	n <- length(a[a$group == "upregulated","group"])
	a$group[a$group=="upregulated"] <- paste("upregulated", " (n=",n,")",sep="")
	#n <- length(a[a$group == "-","group"])
	#a$group[a$group=="-"] <- paste("-unregulated", " (n=",n,")",sep="")

}
#lbls <- as.vector(unlist(a[lblsn]))
#lbls[lbls!=labelName] <- ""

if (is.null(labelName)) {
	a$lbl <- F;
}	else {
	a$lbl <- as.vector(unlist(a[lblsn])) == labelName
}
is <- grep(cargs[2],names(a))
names(a) <- gsub(cargs[2],"",names(a))
head(a,2)

clrs <- c("black","red","steelblue1")

plot.new();pdf(paste(cargs[1],"_splom.pdf",sep=''),h=8,w=8)
#plot.new();png(paste(cargs[1],"_splom.png",sep=''),h=720,w=720)
article.theme$superpose.symbol <- list(col=clrs)
article.theme$superpose.line   <- list(col=clrs,alpha=article.alpha)
print(unique(a$group))
head(a[is])
#splom(~(a[is]),data=a, 
splom(~log2(a[is]),data=a, 
	type = c("g","p", "smooth"),groups=group,lbls=a$lbl,
	par.settings = article.theme,cex=article.cex,pch=".",lwd=article.lwd,font=article.font,fontfamily=article.fontfamily,
	auto.key=list(columns=1,space="top",points=F,col=clrs,font=article.font,fontfamily=article.fontfamily),
	upper.panel = function (x,y,lbls,lwd,...) {
		panel.abline(a=-2, b=1, col="black",lwd=lwd,lty="dashed")
		panel.abline(a=2, b=1, col="black",lwd=lwd, lty="dashed")
		panel.splom(x,y,...)
		panel.xyplot(x[lbls],y[lbls],pch=19,col="red",cex=1,alpha=0.5)
	},
	lower.panel = function (x,y,...) {
 		panel.key(c(sprintf('N=%d\nR=%.2f', length(x),cor(x,y, method="spearman") )),...,col='black',points=F,corner=c(0.5,0.5))
	},
)
dev.off()

