library(lattice)
source(file.path("/import/bc2/home/zavolan/bilebi00/_ARE/scripts/","theme.R"))

cargs <- commandArgs(trailingOnly=T)
#cargs <- c("bins.txt")
a <- read.table(cargs[1], header=F);
names(a) <- c("type","Position","CrosslinkCount","Percent","CumulativePercent", "Length")
head(a)

plot.new(); pdf(paste(cargs[1],".pdf",sep=""),h=8,w=8);
xyplot(CrosslinkCount ~ Position | as.factor(Length), groups=type, data=a, 
	xlab="Position",
	ylab="Crosslink Count",
	main=cargs[1],
	auto.key=list(corner=c(1,1), cex=1.5,font=2),
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=1.5,font=2),
	par.settings = poster.theme,
	lwd=2,cex=1.5, font=2,
	type=c("g","p","l"), span=0.1,
	layout=c(1,length(unique(a$type))))
dev.off()

