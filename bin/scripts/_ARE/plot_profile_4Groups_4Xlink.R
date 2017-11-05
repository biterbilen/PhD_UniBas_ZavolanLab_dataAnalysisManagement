cargs <- commandArgs(trailingOnly = T);
#cargs <- c("slopl500-slopr500.crosslink_distance")

library(latticeExtra)
a <- read.table(cargs[1], header=F);
names(a) <- c("Distance", "Name", "PosteriorSum","Count");
head(a)

a$SiteCount <- sub("\\w+.","",a$Name, perl=T);
a$proteinName <- sub(".site.\\d.\\d","",a$Name, perl=T);
a[grep("hnRNPQ$", a$proteinName, perl=T),"proteinName"] <- "hnRNPQ1";

b <- a[grep("site.[123].", a$SiteCount),];
a <- b;
dim(b)

source(file.path("/import/bc2/home/zavolan/bilebi00/_ARE/scripts/","theme.R"))

plot.new(); pdf(paste(cargs[1], '.PosteriorAverage.SiteCount.pdf', sep=''), h=8, w=8);
xyplot(PosteriorSum/Count ~ Distance | SiteCount, groups=proteinName, data=a, 
	index.cond=list(c(4,5,6,1,2,3)), #this changes the order of the panels
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=1.5,font=2),
	par.settings = poster.theme,
	ylab = "Average Crosslink Score",
	xlab = "Distance from PolyA Site",
	type=c("g","smooth"), span=0.02, 
	auto.key=list(columns=1, lines=T, points=F, corner=c(1,0.90),cex=1.5, lwd=2, font=2)) +
	layer(panel.abline(v=0, col="darkgray",lwd=2, lty=2)) +
	layer(panel.abline(a=0, b=0, col="black",lwd=2))
dev.off()

plot.new(); pdf(paste(cargs[1], '.PosteriorAverage.proteinName.pdf', sep=''), h=8, w=8);
xyplot(PosteriorSum/Count ~ Distance | proteinName, groups=SiteCount, data=a, 
	index.cond=list(c(3,5,1,2,4,6)), #this changes the order of the panels
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=1.5,font=2),
	par.settings = poster.theme,
	ylab = "Average Crosslink Score",
	xlab = "Distance from PolyA Site",
	type=c("g","smooth"), span=0.02, 
	auto.key=list(columns=1, lines=T, points=F, corner=c(1,0.90),cex=1.5, lwd=2, font=2)) +
	layer(panel.abline(v=0, col="darkgray",lwd=2, lty=2)) +
	layer(panel.abline(a=0, b=0, col="black",lwd=2))

dev.off()

quit();
