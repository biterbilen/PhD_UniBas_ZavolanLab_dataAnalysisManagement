library(lattice)
source(file.path("/import/bc2/home/zavolan/bilebi00/_ARE/scripts/","theme.R"))

cargs <- commandArgs(trailingOnly=T);

#cargs <- c("closest.distance.all")

a <- read.table(cargs[1],header=F,sep="\t");

names(a)
names(a) <- c("distance", "frequency", "freq_mean", "freq_sd", "freq_count");
head(a)

b <- NULL;
b$distance <- a$distance;
b$frequency <- a$frequency;
b$type <- "original";
b <- data.frame(b)
head(b)

d <- NULL;
d$distance <- a$distance;
d$frequency <- a$freq_mean + 3*a$freq_sd;
d$type <- "upper_bound";
d <- data.frame(d)

e <- NULL;
e$distance <- a$distance;
e$frequency <- a$freq_mean - 3*a$freq_sd;
e$type <- "lower_bound";
e <- data.frame(e)

k <- rbind(b,d,e)
head(k)

plot.new(); pdf(paste(cargs[1],'.pdf',sep=''));
xyplot(frequency ~ distance, data=k, subset=distance<100, groups=type, 
	ylab="Crosslinked Positions",
	xlab="Distance",
	main=cargs[1],
	auto.key=list(corner=c(1,1), cex=1.5,font=2),
	strip = strip.custom(style=1,bg="gray"),
	par.strip.text = list(cex=1.5,font=2),
	par.settings = poster.theme,
	lwd=2,cex=1.5, font=2,
	type=c("g","l"), span=0.1);
	#type=c("g","smooth"), span=0.1);
dev.off()
