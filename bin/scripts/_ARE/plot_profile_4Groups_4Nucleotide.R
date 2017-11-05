cargs <- commandArgs(trailingOnly = T);


#cargs <- c("nucleotide_distance")
a <- read.table(cargs[1],header=F);
names(a) <- c("Distance", "Name", "Count");

a$SiteCount <- sub("\\w+.","",a$Name, perl=T);
a$Nucleotide <- sub(".site.\\d.\\d","",a$Name, perl=T);

b <- a[grep("site.[123]", a$SiteCount),];

d <- aggregate(Count ~ Distance + SiteCount, sum, data=b)

a <- merge(b, d, by=c("Distance","SiteCount"), suffixes=c("",".all"))
head(a)

library(latticeExtra)

source(file.path("/import/bc2/home/zavolan/bilebi00/_ARE/scripts/","theme.R"))

plot.new(); pdf(paste(cargs[1], '.SiteCount.pdf', sep=''), h=8, w=8);
xyplot(100 * Count / Count.all ~ Distance | SiteCount, groups=Nucleotide, data=a, 
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=1.5,font=2),
	par.settings = poster.theme,
	ylab = "Nucleotide Composition (%)",
	xlab = "Distance from PolyA Site",
	type=c("g","smooth"), span=0.02, 
	auto.key=list(columns=1, lines=T, points=F, corner=c(1,0.90), cex=1.5, lwd=2, font=2)) +
	layer(panel.abline(v=0, col="darkgray",lwd=2, lty=2)) +
	layer(panel.abline(a=25, b=0, col="black",lwd=2))
dev.off()

plot.new(); pdf(paste(cargs[1], '.Nucleotide.pdf', sep=''), h=8, w=8);
xyplot(100 * Count / Count.all ~ Distance | Nucleotide, groups=SiteCount, data=a, 
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=1.5,font=2),
	par.settings = poster.theme,
	ylab = "Nucleotide Composition (%)",
	xlab = "Distance from PolyA Site",
	type=c("g","smooth"), span=0.02, 
	auto.key=list(columns=1, lines=T, points=F, corner=c(1,0.90), cex=1.5, lwd=2, font=2)) +
	layer(panel.abline(v=0, col="darkgray",lwd=2, lty=2)) +
	layer(panel.abline(a=25, b=0, col="black",lwd=2))

dev.off()




quit()
