library(lattice);
source(file.path("/import/bc2/home/zavolan/bilebi00/_ARE/scripts/","theme.R"))


cargs <- commandArgs();

cargs <- c(1,1,1,"all.reg_cor")
a <- read.table(cargs[4], header=F);
names(a) <- c("lib", "d", "count");
dim(a)

libs <- length(unique(a$lib))
if (libs < 3) {
	cols <- libs
} else {
	cols <- 2
}
rows <- (libs %/% cols) + (libs %% cols)

plot.new(); pdf(paste(cargs[4], ".pdf", sep=""), h=8, w=8);
xyplot(count ~ d | as.factor(lib), data=a, 
	xlab="Distance", ylab="Count",
	par.settings = poster.theme,
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=1.5,font=2),
	lwd=2,cex=1.5, font=2,
	layout=c(cols,rows),
	plot.points=F, type=c("g","h"))
dev.off();


