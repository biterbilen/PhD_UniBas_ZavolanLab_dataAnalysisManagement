cargs <- commandArgs(trailingOnly = T);
#cargs <- c("~bilebi00/_EWSR1/Analysis/RNAseq_exonIntronCoverage/all.erfc.introns", "erfc.introns", "siCTRL_1","siEWSR1_1","siCTRL_2","siEWSR1_2","undiff_osteo_P5","20daydiff_osteo_P5");

a <- read.table(cargs[1], header=T,sep="\t")
b <- a[grep(cargs[2], names(a))];
names(b) <- cargs[3:length(cargs)]
print(names(b))
d <- reshape(b, times=names(b), timevar="type", varying=list(names(b)), direction="long",v.names=cargs[2])
head(d)

library(latticeExtra)
library(hexbin)


plot.new(); pdf(paste(cargs[1],'.ecdf.pdf', sep=''), h=8, w=8)
ecdfplot( ~erfc.introns, data=d, groups=type, plot.points=F, 
	  par.strip.text = list(cex=1.5,font=2),
		par.settings=list(axis.text=list(font=2,cex=1.5),
		par.ylab.text=list(font=2,cex=1.5),
		par.xlab.text=list(font=2,cex=1.5),
		par.main.text=list(font=2,cex=1.5)),
	lwd=2,cex=1,font=2,
	na.rm=T, auto.key=list(columns=1, corner=c(0,0.5)), main=basename(cargs[1]))
dev.off();


plot.new(); pdf(paste(cargs[1],'.density.pdf', sep=''), h=8, w=8)
densityplot( ~erfc.introns, data=d, groups=type, plot.points=F, 
	  par.strip.text = list(cex=1.5,font=2),
		par.settings=list(axis.text=list(font=2,cex=1.5),
		par.ylab.text=list(font=2,cex=1.5),
		par.xlab.text=list(font=2,cex=1.5),
		par.main.text=list(font=2,cex=1.5)),
	lwd=2,cex=1,font=2,
	na.rm=T, auto.key=list(columns=1, corner=c(0,0.5)), main=basename(cargs[1]))
dev.off();


plot.new(); pdf(paste(cargs[1], '.scatter.pdf',sep=''), h=8, w=8)
splom(b, main=cargs[1], varname.cex=.8, axis.text.cex=.8, type=c("g","p"),
	par.strip.text = list(cex=1.5,font=2),
	pair.settings=list(axis.text=list(font=2,cex=1.5),
		par.ylab.text=list(font=2,cex=1.5),
		par.xlab.text=list(font=2,cex=1.5),
		par.main.text=list(font=2,cex=1.5)),
	lwd=2,cex=1,font=2,
	lower.panel = function(x, y, ...){
		panel.hexbinplot(x, y, ...)#,
		panel.abline(a=0,b=1, ..., col = 'darkgray')#,
		panel.loess(x, y, ..., col = 'black')#,
		panel.key(c(sprintf('N=%d\nR=%.2f', length(x),cor(x,y, method="spearman") )),...,col='black')
	}, 
	)
dev.off()

alt <- "two.sided"
s1s <- c("20daydiff_osteo_P5","siEWSR1_1", "siEWSR1_2", "siEWSR1_1", "siCTRL_1");
s2s <- c("undiff_osteo_P5"   ,"siCTRL_1",  "siCTRL_2",  "siEWSR1_2", "siCTRL_2"); 

for (i in 1:length(s1s)) {
	s1 <- s1s[i];
	s2 <- s2s[i]
	alt <- "less"
	tt <- t.test(b[c(s1)],b[c(s2)],alternative=alt)
	write.table(format(digits=2, data.frame(list(method="Two sample t-test for erfc.introns ", alternative=alt, sample1=s1, sample2=s2, pvalue=tt$p.value, statistic=tt$statistic))),sep="\t",row.names=F, quote=F);
	alt <- "greater"
	tt <- t.test(b[c(s1)],b[c(s2)],alternative=alt)
	write.table(format(digits=2, data.frame(list(method="Two sample t-test for erfc.introns ", alternative=alt, sample1=s1, sample2=s2, pvalue=tt$p.value, statistic=tt$statistic))),sep="\t",row.names=F, quote=F);
}

#t.test(1:100,y=c(7:200)) 
quit()

