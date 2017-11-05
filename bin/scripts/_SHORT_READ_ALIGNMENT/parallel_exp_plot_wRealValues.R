cargs <- commandArgs();
library(lattice)

#cargs <- c(1,1,1,"../Stepanka/Parallel_plots/dima_clip_targets_in_Hela.clipzQNgene");
#cargs <- c(1,1,1,"../Stepanka/Parallel_plots/DIS3s.clipzQNgene");
b <- read.table(cargs[4], header=T);
a <- b[,-1];

#getpc <- function(x,cs) { 
#	colmin <- 1:cs * 0;
#	for (c in 1:cs) {
#		colmin[c] <- sort(x[x[,c]>0,c])[1] #get min that is greater than 0
#	};
#	return(colmin);
#}; 

rs <- dim(a)[1]
cs <- dim(a)[2]-1 
#mins <- getpc(a,cs);
#d <- a[,1:cs] + matrix(rep(mins,rs),nrow=rs,byrow=T)
d <- a[,1:cs] + 1;

a[,1:cs] <- log2(d)

rng <- range(unlist(lapply(lapply(a, as.numeric), range)))
prng <- pretty(rng)

ltheme <- canonical.theme(color = FALSE)      ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme, lwd=4)      ## set as default

plot.new(); pdf(paste(cargs[4],'.realExp.parallel.pdf', sep=''),h=8,w=8);
parallel(~a[,1:cs]|gid, data=a, xlab="Log2 Expression",common.scale=T, main=cargs[4],
	par.settings=list(axis.text=list(font=2,cex=1.0),
		par.strip.text=list(font=2,cex=1.0),
		par.ylab.text=list(font=2,cex=1.0),
		par.xlab.text=list(font=2,cex=1.0),
		par.main.text=list(font=2,cex=1.0)),
	lwd=4,
	type=c("g","l"),
	layout=c(1,3),
	scales = list(x = list(rot=90, at = (prng-min(rng))/diff(rng), labels = prng)))
dev.off();

#plot.new(); pdf(paste(cargs[4],'.parallel.pdf', sep=''),h=8,w=8);
#parallel(~ asinh(a[c(2:12)]), data = a, alpha = 0.1, lty = 1) 
#dev.off()

