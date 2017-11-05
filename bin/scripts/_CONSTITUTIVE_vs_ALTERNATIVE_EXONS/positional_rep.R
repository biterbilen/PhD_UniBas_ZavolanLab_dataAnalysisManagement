cargs <- commandArgs();
print (args)

library("latticeExtra") 

#cargs=c(1,1,1,"c");
a <- read.table(cargs[4], header=T, check.names=F)

tit = paste("Positional Reproducibility of Xlink Sites", sep="")
ofilename <- paste(cargs[4], ".pdf", sep="")
plot.new(); pdf(ofilename, width=8, height=8)

nop <- xyplot(log(numberOfPositions) ~ (minT2C) | proteinName,
	data = a, type = c("smooth", "g"), groups=sampleId,
	#auto.key=list(title=tit, columns=5, layout=c(1,10)))
	auto.key=list(title=tit))
forp <- xyplot(fracOfReproducedPositions ~ (minT2C)| proteinName, 
	data = a, type = c("o"), groups=sampleId)
doubleYScale(forp, nop, style1 = 0, style2 = 3, add.ylab2 = T)

#xyplot(log(numberOfPositions)+fracOfReproducedPositions ~ log(minT2C) | proteinName, data=a, type="l", groups=sampleId, auto.key=list(title=tit, columns=3), layout=c(2,5)) 

dev.off()
