library(latticeExtra);

#TES and TSS
cargs <- c("ENCODEHekRep3.INTRON4TR", "1","INTRON")
cargs <- c("ENCODEHekRep3.TES4TR", "1","TES")
cargs <- commandArgs(trailingOnly=T)
minval <- 1
if (length(cargs) > 1) {
	minval <- as.numeric(cargs[2])
}
panel <- "panel"
if (length(cargs) > 2) {
	panel <- cargs[3]
}
print(minval)
print(panel)

#pausing index
b <- read.table(cargs[1],header=T,stringsAsFactors=F);
head(b,2)
b$panel <- panel
a <- b[b$region.POL2Chip>minval & b$gene.POL2Chip>minval,]
dim(b)
dim(a)
head(a,3) 
a$group <- "other";
a[a$X0>0  & a$X0.1==0 & a$X0.2==0,]$group <- "EWSR1" 
a[a$X0==0 & a$X0.1>0  & a$X0.2==0,]$group <- "Rloop" 
a[a$X0==0 & a$X0.1==0 & a$X0.2>0, ]$group <- "Cstf64" 
a[a$X0 + a$X0.2 == 2, ]$group <- "EWSR1.Cstf64" 
a[a$X0 + a$X0.1 == 2, ]$group <- "EWSR1.Rloop" 
#a[a$X0 + a$X0.1 + a$X0.2==3, ]$group <- "EWSR1.Rloop.Cstf64" 
#k <- a[a$group == "both",]; k$group <- "EWSR1"
#a[a$group == "both",]$group <- "Rloop";
#a <- rbind(k,a);
xtabs(~a$group)

a$region.POL2.TR <- (a$region.POL2Chip/a$gene.POL2Chip) * (a$gene.length/a$region.length)
d <- a[c("panel","group","region.POL2.TR","id","gid")];
d$panel <- paste(panel,sep="")
write.table(d,file=paste(cargs[1],".out",sep=""),sep="\t",quote=F,row.names=F)

#ecdfplot(~ TES.POL2.TR, group=group, data=a, 
#	scales=list(x=list(log=F),y=list(log=F)),
#	types=c("l","g"),
##	subset=length>1000,
#	#subset=aroundINTRON.length>300&INTRON.length>10000,
#	auto.key=list(space="top"))

#-------------------------------------------------------------------

