library(latticeExtra)

cargs <- c("distance_count")
cargs <- commandArgs(trailingOnly=T)

a <- read.table(cargs[1],header=F)
names(a) <- c("distance","count")
plot.new(); pdf(paste(cargs[1],".pdf",sep=""),h=8,w=8)
xyplot(cumsum(a$count)/sum(a$count) ~ distance,data=a,type=c("g","l"),scales=list(x=list(log=T),x=list(log=T)))
dev.off()
