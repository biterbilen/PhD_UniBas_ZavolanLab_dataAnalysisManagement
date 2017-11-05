library(latticeExtra)

cargs <- commandArgs(trailingOnly=T);
#cargs <- c("RNAseq.raw.qn.w_geneSymbol.FC.w_Hancock.w_Gaussian.w_outlier.w_Kishore")

a <- read.table(cargs[1], header=T);
head(a)


plot.new(); pdf(paste(cargs[1],".ellipse_groups.pdf",sep=''),h=8,w=8);
xyplot(mRNAseq_EWSR1_siEWSR1_1_over_mRNAseq_EWSR1_siCTRL_1 ~ mRNAseq_EWSR1_siEWSR1_2_over_mRNAseq_EWSR1_siCTRL_2 | is_outlier, 
	layout=c(1,2),
	groups=Hancock2008Regulation, 
	type=c("g","p"),
	data = a, scales = "free",
	par.settings = list(superpose.symbol = list(pch=c(15:17),cex=c(0.5,0.15,0.5),
		superpose.line = list(lwd=c(4,4,4), lty=2:4))),
	panel = function(x, y, ...) {
		panel.xyplot(x, y, ...)
#		panel.ellipse(x, y, robust=T,...)
	},
	auto.key = list(space="top",columns=3))
dev.off()

quit();

plot.new(); pdf(paste(cargs[1],".ellipse_all.pdf",sep=''),h=8,w=8);
xyplot(mRNAseq_EWSR1_siEWSR1_1_over_mRNAseq_EWSR1_siCTRL_1 ~ mRNAseq_EWSR1_siEWSR1_2_over_mRNAseq_EWSR1_siCTRL_2, 
	data = a, scales = "free",
	par.settings = list(superpose.symbol = list(pch=c(15:17),cex=c(0.35,0.1,0.35),
		superpose.line = list(lwd=c(4,4,4), lty=2:4))),
	panel = function(x, y, ...) {
		panel.xyplot(x, y, ...)
		panel.ellipse(x, y, robust=F,col="red",...)
		#panel.ellipse(x, y, robust=T,col="red",...)
	},
	auto.key = list(x = .1, y = .8, corner = c(0, 0)))
dev.off()

library(MASS)
d <- a[grep("mRNA",names(a))]
dd <- cov.trob(d)
head(dd)
dd$n.obs
dd$cov
dd$wt
dd$center

x <- d[,1]
y <- d[,2]
kde2d(x,y, h, n = 25, lims = c(range(x), range(y)))
