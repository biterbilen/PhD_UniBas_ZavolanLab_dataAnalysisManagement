cargs <- commandArgs();
cargs=c(1,1,1,"../PAS_count");
a <- read.table(cargs[4], header=T);

library(lattice)
cor.a <- cor(as.matrix(a),use="pair")
pl <- splom(a)
print(pl)

library(ggplot2)
plot.new(); pdf(paste(cargs[4],".pdf",sep=""), w=8,h=8);
plotmatrix(a)
#plotmatrix(a, auto.key=list(title="PAS per gene"))
dev.off();


