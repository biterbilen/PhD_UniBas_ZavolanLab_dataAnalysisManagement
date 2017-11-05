library(lattice)
cargs <- commandArgs()
#cargs <- c(1,1,1,"a", "DIS3L2");
ext <- '.cor';
a <- read.table(cargs[4], header=T)
rs <- dim(a)[1];
cs <- dim(a)[2];
index <- which(a$gid == cargs[5]);
D <- log2(a[,3:cs-1]+1)
dis3l2 <- t(D[index,])
corr <- apply(D, 1, function(x) { cor(dis3l2,x); })
corr[index]
ocorr <- data.frame(a$id, a$gid, corr);
names(ocorr) <- c("id", "gid", "cor");
socorr <- ocorr[order(corr),];
outfile=paste(cargs[4], ext, sep="");
write.table(socorr, file=outfile, sep="\t", row.names=F, quote=F)

plot.new(); pdf(paste(cargs[4], ext, '.cor.pdf', sep=''), h=8, w=8);
#plot(1:rs-2,socorr$cor)
densityplot(socorr$cor, type=c("percent"));
dev.off()
  
