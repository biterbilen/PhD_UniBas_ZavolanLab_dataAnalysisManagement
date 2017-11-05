library(lattice)

cargs <- commandArgs(trailingOnly=T)
#cargs <- c("RegionExpression.CTRL_KD.feature.stats",1)

#a <- read.table("features_significance_intronHEK293/cons",header=T)
a <- read.table(cargs[1],header=T,sep="\t")
head(a)
paired <- F;
alt <- "greater"
if (length(cargs) > 1) {
	paired <- as.logical(as.numeric(cargs[2]));
}
if (length(cargs) > 2) {
	alt <- cargs[3];
}

libs <- unique(a$lib)
tags <- sort(unique(a$tag))
types <- unique(a$type)
paste(libs)
paste(tags)
paste(types)

for (type in types) {
	print (type)
	for (lib in libs) {
		print(lib)
		x <- as.numeric(a[a$type == type & a$lib == lib & a$tag == tags[1],"value"]);
		y <- as.numeric(a[a$type == type & a$lib == lib & a$tag == tags[2],"value"]);
		print(c(type, lib, length(x), length(y)))
		tt  <- t.test(x,y,alt="greater",paired=paired)

		write.table(
			format(digits=2, data.frame(list(method=paste("Two sample t-test (paired=", as.numeric(paired), ")", sep=""), type=type, lib=lib, alternative=alt, sample1=tags[1], sample2=tags[2], pvalue=tt$p.value, statistic=tt$statistic,sep="----------", sample1mean=mean(x),sample2mean=mean(y),sample1stdev=sd(x),sample2stdev=sd(y),sample1count=length(x),sample2count=length(y)))),
			sep="\t",row.names=F, quote=F);
	}
}

