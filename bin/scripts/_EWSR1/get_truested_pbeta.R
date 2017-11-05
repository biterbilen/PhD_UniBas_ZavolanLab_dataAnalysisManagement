#cargs <- c(1,1,1,"_288_289.introns.xlinkT_mutrate0.08.txt");
cargs <- commandArgs()

a <- read.table(cargs[4],header=T);

a[is.na(a$pbeta), c("pbeta")] <- 1 #XXX better way of assigning pbeta nonexisting data
a[a$pbeta==0, c("pbeta")] <- .Machine$double.xmin

a[is.na(a$pbeta.1), c("pbeta.1")] <- 1 #XXX better way of assigning pbeta nonexisting data
a[a$pbeta.1==0, c("pbeta.1")] <- .Machine$double.xmin

#--
head(sort(a[is.na(a$freq), c("pbeta.1")],decreasing=F),2)
head(sort(a[is.na(a$freq.1), c("pbeta")], decreasing=F),2)

#minimum pbetas that is not reproduced in between the libraries
#get best score to be in background
minpbeta <- min(c(min(a[is.na(a$freq), c("pbeta.1")]), min(a[is.na(a$freq.1), c("pbeta")])))
print(minpbeta)

#get the maximum pbeta as a background score; ie min score to be in the foreground
a$maxpbeta <- apply(a[c("pbeta","pbeta.1")], 1, max)

#trusted set is the ones that is better then the best score to be in the background
b <- (a[a$maxpbeta < minpbeta,])

#--
d1 <- sort(b$maxpbeta, decreasing=F,index.return=T);
d2 <- sort(a$maxpbeta, decreasing=F,index.return=T);
fields <- c("seqname", "start", "end", "strand", "gene_id.x", "maxpbeta");
write.table(b[d1$ix,fields], col.names=gsub(".x", "", fields, fixed=T), row.names=F, quote=F, sep="\t", file=paste("trusted", cargs[4],sep=''))
write.table(a[d2$ix,fields], col.names=gsub(".x", "", fields, fixed=T), row.names=F, quote=F, sep="\t", file=paste("full", cargs[4],sep=''))

#--
library(lattice)
library(hexbin)

xs <- 2;
ys <- 1;
p1 <- xyplot(log2(freq)  ~ log2(freq.1), data=a, panel=panel.hexbinplot);
p2 <- xyplot(log2(pbeta) ~ log2(pbeta.1), data=a, panel=panel.hexbinplot, subset=pbeta>0 & pbeta.1>0);
plot.new(); pdf(paste(cargs[4], "_params.pdf", sep=''), h=8, w=8);
plot(p1, split = c(1, 1, xs, ys))
plot(p2, split = c(2, 1, xs, ys), newpage = FALSE)
dev.off();

