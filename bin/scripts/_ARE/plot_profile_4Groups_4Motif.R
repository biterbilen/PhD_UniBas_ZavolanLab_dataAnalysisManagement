cargs <- commandArgs(trailingOnly = T);

#cargs <- c("TIA1.motif_distance","3","2");
rows <- 2;
cols <- 3;
if (length(cargs) > 1) {
	cols <- as.numeric(cargs[2]);
}
if (length(cargs) > 2) {
	rows <- as.numeric(cargs[3]);
}
print(cols)

a <- read.table(cargs[1],header=F,sep="\t");
names(a) <- c("Distance", "proteinName", "Count");

d <- aggregate(Count ~ proteinName, sum, data=a)

b <- merge(a, d, by=c("proteinName"), suffixes=c("",".all"))
head(b)

library(latticeExtra)

source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))

plot.new(); pdf(paste(cargs[1], '.crosslink.pdf', sep=''), h=8, w=8);
xyplot(100 * Count / Count.all ~ Distance | as.factor(proteinName), data=b, 
	layout=c(rows,cols),
	scales=list(x=list(relation="free")),
	strip = strip.custom(style=1,bg="white"),
	par.strip.text = list(cex=article.cex,font=article.font),
	par.settings = article.theme,
	ylab = "Motif Start (%)",
	xlab = "Distance to Crosslink Site",
	type=c("g","l","p"),
	col = "black", lwd=article.lwd, 
#	auto.key=list(columns=1, lines=T, points=F, corner=c(1,0.90), cex=1.5, lwd=2, font=2)
	) +
#	layer(panel.abline(a=0.05, b=0, col="black",lwd=2)) +
	layer(panel.abline(v=0, col="darkgray",lwd=2, lty=2))

dev.off()

quit()
