library(lattice)
source(file.path("/import/bc2/home/zavolan/bilebi00/_CLIP/scripts/","theme.R"))

cargs <- commandArgs(trailingOnly=T);
#cargs <- c("/import/bc2/home/zavolan/bilebi00/_CLIP/Analysis/Project_DIS3L2/xlinkEnrichment/tRNAProfile_L72/72.tRNAs.", "title");
ylab <- "Crosslinked Positions";

if (length(cargs) > 1) {
	ylab <- cargs[2];
}

a <- read.table(cargs[1],header=F,sep="\t");

names(a)
names(a) <- c("distance", "frequency", "freq_min", "freq_max", "freq_count");
head(a)

b <- NULL;
b$distance <- a$distance;
b$frequency <- a$frequency;
b$type <- "original";
b <- data.frame(b)
head(b)

d <- NULL;
d$distance <- a$distance;
d$frequency <- a$freq_max
d$type <- "upper CI";
d <- data.frame(d)

e <- NULL;
e$distance <- a$distance;
e$frequency <- a$freq_min
e$type <- "lower CI";
e <- data.frame(e)

k <- rbind(b,d,e)
head(k)

plot.new(); pdf(paste(cargs[1],'.pdf',sep=''));
xyplot(frequency ~ distance, data=k, subset=distance<100, groups=type, 
	ylab=ylab,
	xlab="Distance",
	main=cargs[1],
	auto.key=list(space="top", columns=3,cex=article.cex,font=article.font,points=F,lines=F,col=c(article.colors[1])),
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=article.cex,font=article.font),
	par.settings = article.theme, col=c(article.colors[1]),
	lty=c("solid","dashed","dashed"),
	lwd=c(article.lwd,article.lwd/2,article.lwd/2),cex=article.cex, font=article.font,
	type=c("g","l"), pch=16);
	#type=c("g","smooth"), span=0.1);
dev.off()
