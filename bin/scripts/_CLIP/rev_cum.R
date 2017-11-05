library(latticeExtra)
source(file.path("/import/bc2/home/zavolan/bilebi00/_ARE/scripts/","theme.R"))

cargs <- commandArgs(trailingOnly=T)
#cargs <- c("out","-xlinkEnrichment","1", "2","a","b","EWSR1 288","EWSR1 289")

filecount <- as.numeric(cargs[4])
d <- NULL;
for (i in 1:filecount) {
	print(i)
	print(cargs[4+i])
	print(4+i+filecount)
	print(cargs[4+i+filecount])
	a <- read.table(cargs[4+i]);
	a$type  <- cargs[4+i+filecount]
	d <- rbind(d,a)
}
if (cargs[3] == "1") {
	d$V1 = -d$V1;
}
head(d)

ylabel <- "Reverse Cumulative"
if (length(cargs) == 4 + (2*filecount) + 1) {
	ylabel <- cargs[4 + (2*filecount) + 1];
}	

print(ylabel)

plot.new(); pdf(paste(cargs[1],'.pdf',sep=''),h=8,w=8);
xyplot(V2 ~ V1, groups=type, data=d, type=c("g","p"),
	xlab = cargs[2],ylab=ylabel,
	par.settings = article.theme,
	xscale.components = xscale.components.log10ticks,
	yscale.components = yscale.components.log10ticks, axis = axis.grid,
	scales=list(x=list(log=T,tck=c(1,1)),y=list(log=T)), auto.key=list(space="top",cex=article.cex))
dev.off();

