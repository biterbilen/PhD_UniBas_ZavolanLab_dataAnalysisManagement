cargs <- commandArgs();
#TODO write
library(lattice)
cargs <- c(1,1,1,"/import/bc2/home/zavolan/bilebi00/_DIS3/Normalization_PAPD5_25nc/","prenorm_reverse_cum", "^raw*");

files <- dir(cargs[4], pattern=cargs[6]);

plot.new(); pdf(paste(cargs[4],cargs[5],'.pdf',sep=""), w=8, h=8);
i <- 1;
plots= c();
for (file in files) {
	a <- read.table(paste(cargs[4],file,sep=''), header=F);
	p <- xyplot(log2(V2) ~ log2(V1), data=a, xlab="log2(t)", ylab="log2(# of (t > T))", style=i)
	if ( i == 1) { 
		plot(p);
	else {
		plot(p, newpage=F);
	}
	plots <- cbind(plots,c(p));
	i <- i + 1;
}

update(plots)
dev.off();

