cargs <- commandArgs();

library(latticeExtra)

cargs <- c(1,1,1,"Normalization_PAPD5/table.norm_and_raw")
#cargs <- c(1,1,1,"Normalization_PAPD5_T2C/table.norm_and_raw")
aa <- read.table(cargs[4],header=T);

copies="copies"
copies=""
xi=5
yi=8
r=5
r=0
outf=paste(cargs[4],".raw",copies,".PAPD5_2-IGF2BP1.pdf",sep="");
outf=paste(cargs[4],".norm",copies,".PAPD5_2-IGF2BP1.pdf",sep="");

xi=4
yi=8
r=5
r=0
outf=paste(cargs[4],".raw",copies,".PAPD5_1-IGF2BP1.pdf",sep="");
outf=paste(cargs[4],".norm",copies,".PAPD5_1-IGF2BP1.pdf",sep="");

xi=4
yi=5
r=5
r=0
outf=paste(cargs[4],".raw",copies,".PAPD5.pdf",sep="");
outf=paste(cargs[4],".norm",copies,".PAPD5.pdf",sep="");

xi=6
yi=7
r=5
r=0
outf=paste(cargs[4],".raw",copies,".TIA1.pdf",sep="");
outf=paste(cargs[4],".norm",copies,".TIA1.pdf",sep="");

sel=apply(as.matrix(cbind(aa[,xi+r], aa[,yi+r])), 1, function (x) { x[1] > 0 && x[2] > 0 });
a <- aa[sel,]

outf
plot.new(); pdf(outf,h=8,w=8); 
xyplot(log2(a[,yi+r]) ~ log2(a[,xi+r]), data=a, main=outf) +
 layer(panel.key(c(sprintf('R=%.2f N=%d', cor(x,y), length(x))),lwd=2,cex=2));
dev.off()



