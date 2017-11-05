source(file.path("/import/bc2/home/zavolan/bilebi00/_EWSR1/scripts/","getDataFromGtfPlus.R"))
library(lattice)

cargs <- commandArgs(trailingOnly=T)
#cargs <- c("RegionExpression.CTRL_KD.feature.stats",1)

#a <- read.table("features_significance_intronHEK293/cons",header=T)
#a <- read.table(cargs[1],header=T,sep="\t")
a <- getDataFromGtfPlus(cargs[1], c("effective_length","value"))
head(a,3)
a$lib <- sub ("\\..+", "",a$source)
a$tag <- sub ("([^.]+).", "",a$source)

sub ("\\..+", "",unique(a$feature))
sub ("([^.]+).", "",unique(a$feature))

head(a, 3)
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
features <- unique(a$feature)
paste(libs)
paste(tags)
paste(features)

for (feature in features) {
	print (feature)
	for (lib in libs) {
		print(lib)
		x <- as.numeric(a[a$feature == feature & a$lib == lib & a$tag == tags[1],"value"]);
		y <- as.numeric(a[a$feature == feature & a$lib == lib & a$tag == tags[2],"value"]);
		print(c(feature, lib, length(x), length(y)))
		tt  <- t.test(x,y,alt="greater",paired=paired)

		write.table(
			format(digits=2, data.frame(list(method=paste("Two sample t-test (paired=", as.numeric(paired), ")", sep=""), feature=feature, lib=lib, alternative=alt, sample1=tags[1], sample2=tags[2], pvalue=tt$p.value, statistic=tt$statistic,sep="----------", sample1mean=mean(x),sample2mean=mean(y),sample1stdev=sd(x),sample2stdev=sd(y),sample1count=length(x),sample2count=length(y)))),
			sep="\t",row.names=F, quote=F);
	}
}

