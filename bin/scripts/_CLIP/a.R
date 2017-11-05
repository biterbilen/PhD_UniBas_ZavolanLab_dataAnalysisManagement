library(latticeExtra)
library(hexbin)

cargs <- c("RNAseq/input_4_bayseq", "RNAseq/lengths")
cargs <- commandArgs(trailingOnly=T)
dir()
a <- read.table(cargs[1],header=T)
len <- read.table(cargs[2],header=F)
len <- len$V1
head(len)
nsamples <- dim(a)[2] - 1
ngenes <- dim(a)[1]
dim(a)
length(len)
head(as.list(len))

counts <- 1000 * a[,-1] / matrix(rep(len, nsamples),byrow=F,nrow=ngenes)
sums <- as.vector(colSums(counts))
fractions <- counts / matrix(rep(sums,ngenes),byrow=T,ncol=nsamples)

write.table(cbind(a[1],fractions), file=paste(cargs[1],'.txt',sep=''), sep="\t",quote=F,row.names=F)

plot.new();pdf(paste(cargs[1],'.pdf',sep=''),h=8,w=8)
splom(log2(fractions),cex=1,pch=".",col="black")
#splom(log2(fractions),cex=2,pch=".",col="black",panel=panel.hexbinplot)
dev.off()

