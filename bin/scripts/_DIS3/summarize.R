cargs = commandArgs();
#TODO change
cargs = c(1,1,1,"~bilebi00/DB/merged_der");
library(lattice) 

a = read.table(cargs[4], header=TRUE, sep="\t");

#annots
annots <- as.data.frame(xtabs(~description+annot_type,a))
tit <- 'Trusted region annotations -DIS3 wrt backgrounds';
ofile <- paste(cargs[4], "_annotations_freq.pdf", sep='');
plot.new(); pdf(ofile, width=8, height=8)
barchart(log2(Freq) ~ annot_type|description, data=annots, stack=T,
	par.strip.text=list(cex=0.5), 
	par.settings=list(axis.text=list(cex=0.3),vertical=T),
	xlab=list(cex=2), ylab=list(cex=2),
	main=list(label=tit, cex=2),
	cex=0.1,
	#lcw=c(0.5),
	auto.key = list(columns = 1, cex=0.5), groups=description,type=c("g"))
dev.off()

#TODO
#this part doesn't do the intended, sum(Freq) is for all groups nor for each group separately
annots <- as.data.frame(xtabs(~description+annot_type,a))
tit <- 'Trusted region annotations -DIS3 wrt backgrounds';
ofile <- paste(cargs[4], "_annotations_frac.pdf", sep='');
plot.new(); pdf(ofile, width=8, height=8)
barchart(Freq/sum(Freq)~ annot_type|description, data=annots, stack=T,
	par.strip.text=list(cex=0.5), 
	par.settings=list(axis.text=list(cex=0.3),vertical=T),
	xlab=list(cex=2), ylab=list(cex=2),
	main=list(label=tit, cex=2),
	cex=0.1,
	#lcw=c(0.5),
	auto.key = list(columns = 1, cex=0.5), groups=description,type=c("g"))
dev.off()

#scatter of zSscore
tit <- 'zScore of DIS3 wrt backgrounds';
ofile <- paste(cargs[4], "_zScoreScatter.pdf", sep='');
plot.new(); pdf(ofile, width=8, height=8)
xyplot(score2 ~ score | description, data=a,
	par.strip.text=list(cex=0.5), 
	#par.settings=list(axis.text=list(cex=1.5)),
	xlab=list(cex=2), ylab=list(cex=2),
	main=list(label=tit, cex=2),
	cex=0.1,
	#lcw=c(0.5),
	auto.key = list(columns = 1, cex=0.5), groups=description,type=c("p","g"))
dev.off()



