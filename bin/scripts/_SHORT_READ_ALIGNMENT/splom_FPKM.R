library(latticeExtra);
library(hexbin);

cargs <- c(1,1,1,"../Scheiffele/Cuffdiff_trusted_RNAseq/genes.fpkm_tracking",0);

a <- read.table(cargs[4], header=T);
cutoff<-cargs[5];
#mygt <- function(x) { x[1]>cutoff && x[2]>cutoff && x[3]>cutoff && x[4]>cutoff && x[5]>cutoff && x[6]>cutoff && x[7]>cutoff && x[8]>cutoff }; 
mygt <- function(x) { x[1]>cutoff && x[2]>cutoff && x[3]>cutoff && x[4]>cutoff };
	#TODO generalize column selection
b <- apply(a[,c(11,14,17,20)],1, mygt); 
d <- a[b,c(11,14,17,20)]; 

#plot.new(); pdf(paste(cargs[4], '.pdf', sep=''), h=8, w=8);
#splom(log(a[2:5]+1), panel=panel.hexbinplot,main=cargs[4],
#plot.new(); pdf(paste(cargs[4], 'cut', cutoff,'.pdf', sep=''), h=8, w=8);
plot.new(); pdf(paste('cut', cutoff,'.pdf', sep=''), h=8, w=8);
splom(log(d), panel=panel.hexbinplot,main=cargs[4], varname.cex=0.35, axis.text.cex=0.5,
	lower.panel = function(x, y, ...){
		panel.hexbinplot(x, y, ...)#,
		panel.loess(x, y, ..., col = 'red')#,
		panel.key(c(sprintf('N=%d; R=%.2f', length(x),cor(x,y) )),...,col='black',cex=0.5)
	},
	)
dev.off();

