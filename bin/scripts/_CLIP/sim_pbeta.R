source(file.path("/import/bc2/home/zavolan/bilebi00/_ARE/scripts/","theme.R"))

library(latticeExtra)
library(hexbin)
library(mclust)

cargs <- commandArgs(trailingOnly = T)

s <- a$SNP == F & a$Contaminated == F & a$Mismapped == F & a$HypoCovered == F & a$HyperMutatedHighlyCovered == F
mu <- mean(a[s,"mut"] / a[s,"copies"])

#Estimate muRef using trusted set and those that are mutated already once
sb1 <- s & a$mut_T2A > 0
muRef1 <- mean(a[sb1,"mut_T2A"] / a[sb1,"copies"])
sb2 <- s & a$mut_T2G > 0
muRef2 <- mean(a[sb2,"mut_T2G"] / a[sb2,"copies"])

#get the maximum
muRef <- max(muRef1, muRef2)

#calculate pbeta
#TODO decide here
aa <- min(a[s,"mut"])
bb <- ceiling(aa / muRef) - aa;
print(write.table(cbind("mus:", mu, muRef1, muRef2), quote=F, row.names=F, col.names=F))
print(write.table(cbind("mulens:", length(s[s==T]), length(sb1[sb1==T]), length(sb2[sb2==T])), quote=F, row.names=F, col.names=F))
print(write.table(cbind("pbeta hyperparameters:", muRef, aa, bb), quote=F, sep="\t", row.names=F, col.names=F))
a$pbeta <- pbeta(muRef, a$mut + aa, a$copies - a$mut + bb, log=T);

#plot scores
select <- sample((1:nrow(a))[s], size=2000)
p <- levelplot(pbeta ~ copies + mut / copies, data=a[select,], 
	col.regions=colorRampPalette(c("red","orange","yellow","gray","black")),
	colorkey = list(space = "top"), main="pbeta Distribution",
	xlab = "Selected Tags", ylab="Position Mutations / Selected Tags",
	cex=article.cex, font=article.font, lwd=article.lwd,
	par.settings = article.theme, type=c("g","p"),
	panel = panel.levelplot.points,
	scales=list(x=list("log" = T)),
	xscale.components = xscale.components.log10ticks)

ofile <- paste(cargs[2], ".pbeta.pdf", sep="")
plot.new(); pdf(ofile, h=8,w=8)
print(p)
dev.off()
