library(latticeExtra);

#TES and TSS
cargs <- c("ENCODEHelas3Rep1.TES4TR", "ENCODEHelas3Rep1.TSS4TR", "0");
cargs <- commandArgs(trailingOnly=T)
minval <- 1
if (length(cargs) > 2) {
	minval <- as.numeric(cargs[3])
}
print(minval)

#TES pausing
b <- read.table(cargs[1],header=T,stringsAsFactors=F);
head(b,2)
b$panel <- "TES.POL2.TR";
a <- b[b$TES.POL2Chip>minval & b$extTES.POL2Chip>minval,]

head(a,2)
a$group <- "other";
#a[a$X3UTR.1+a$X5UTR.1+a$CDS.1+a$DGE.1+a$INTRON.1+a$UGS.1>0,]$group <- "MAYBECLIPPED" 
#a[a$X3UTR+a$X5UTR+a$CDS+a$DGE+a$INTRON+a$UGS>0,]$group <- "CLIPPED" 
#a[a$X3UTR+a$DGE>0,]$group <- "END" 
a[a$X0>0,]$group <- "EWSR1" 
a[a$X0.1>0,]$group <- "Rloop" 
a[a$X0>0&a$X0.1>0,]$group <- "both" 
#k <- a[a$group == "both",]; k$group <- "EWSR1"
#a[a$group == "both",]$group <- "Rloop";
#a <- rbind(k,a);
xtabs(~a$group)

a$TES.POL2.TR <- (a$TES.POL2Chip/a$extTES.POL2Chip) * (a$extTES.length/a$TES.length)
d <- a[c("panel","group","TES.POL2.TR","gid")];
d$panel <- "TravelRateinTES";
write.table(d,file=paste(cargs[1],".out",sep=""),sep="\t",quote=F,row.names=F)

#ecdfplot(~ TES.POL2.TR, group=group, data=a, 
#	scales=list(x=list(log=F),y=list(log=F)),
#	types=c("l","g"),
##	subset=length>1000,
#	#subset=aroundINTRON.length>300&INTRON.length>10000,
#	auto.key=list(space="top"))

#-------------------------------------------------------------------

#TSS pausing
b <- read.table(cargs[2],header=T,stringsAsFactors=F);
b$panel <- "TSS.POL2.TR";
a <- b[b$TSS.POL2Chip>minval & b$extTSS.POL2Chip>minval,]

head(a,2)
a$group <- "other";
#a[a$X3UTR.1+a$X5UTR.1+a$CDS.1+a$DGE.1+a$INTRON.1+a$UGS.1>0,]$group <- "MAYBECLIPPED" 
#a[a$X3UTR+a$X5UTR+a$CDS+a$DGE+a$INTRON+a$UGS>0,]$group <- "CLIPPED" 
#a[a$INTRON>0,]$group <- "INTRON" 
#a[a$X0>0,]$group <- "TSS"
#a[a$X0.1>0,]$group <- "TSS.INTRON"
a[a$X0>0,]$group <- "EWSR1" 
a[a$X0.1>0,]$group <- "Rloop" 
a[a$X0>0&a$X0.1>0,]$group <- "both" 
#k <- a[a$group == "both",]; k$group <- "EWSR1"
#a[a$group == "both",]$group <- "Rloop";
#a <- rbind(k,a);
xtabs(~a$group)

a$TSS.POL2.TR <- (a$TSS.POL2Chip/a$extTSS.POL2Chip) * (a$extTSS.length/a$TSS.length)
d <- a[c("panel","group","TSS.POL2.TR","gid")];
d$panel <- "TravelRateinTSS";
write.table(d,file=paste(cargs[2],".out",sep=""),sep="\t",quote=F,row.names=F)

#ecdfplot(~ TSS.POL2.TR, group=group, data=a, 
#	scales=list(x=list(log=F),y=list(log=F)),
#	types=c("l","g"),
##	subset=length>1000,
#	#subset=aroundINTRON.length>300&INTRON.length>10000,
#	auto.key=list(space="top"))

quit()

#-------------------------------------------------------------------
a <- read.table("all4TR",header=T,stringsAsFactors=F);
b <- read.table("intron4TR",header=T,stringsAsFactors=F);
a <- b[b$KD.INTRON>0 & b$Ctrl.INTRON>0,]
head(a,2)
#a$displacementSe <- a$displacementSd / sqrt(a$displacementCount);
a$KD.TR <- (a$KD.INTRON/a$KD.aroundINTRON) * (a$aroundINTRON.length/a$INTRON.length)
a$Ctrl.TR <- (a$Ctrl.INTRON/a$Ctrl.aroundINTRON) * (a$aroundINTRON.length/a$INTRON.length)
a$KDdivCtrl.TR <- a$KD.TR / a$Ctrl.TR
#a$KD.TR <- (a$KD.INTRON/a$INTRON.length) / (a$KD.aroundINTRON/a$aroundINTRON.length)
#a$Ctrl.TR <- (a$Ctrl.INTRON/a$INTRON.length) / (a$Ctrl.aroundINTRON/a$aroundINTRON.length)
#a$KD.TR <- (a$KD.TSS/400) / (a$KD.TRX/a$length)
#a$Ctrl.TR <- (a$Ctrl.TSS/400) / (a$Ctrl.TRX/a$length)
#a$KD.TR <- a$KD.TSS / a$siEWSR1_1
#a$Ctrl.TR <- a$Ctrl.TSS / a$siEWSR1_1
a$group <- "other" 
a[a$INTRON.xlink>0,]$group <- "MAYBECLIPPED" 
a[a$INTRON.xlink>10,]$group <- "CLIPPED" 
#a[a$X3UTR.1+a$X5UTR.1+a$CDS.1+a$DGE.1+a$INTRON.1+a$UGS.1>0,]$group <- "MAYBECLIPPED" 
#a[a$X3UTR+a$X5UTR+a$CDS+a$DGE+a$INTRON+a$UGS>0,]$group <- "CLIPPED" 
#a[!a$panel=="NULL",]$group <- "FIRST.INTRON" 
e <- a[a$group=="MAYBECLIPPED",c("KDdivCtrl.TR","INTRON.xlink")]
head(e)
splom(log(e))

cor(e$KDdivCtrl.TR,e$INTRON.xlink)
fivenum(e$KDdivCtrl.TR)
fivenum(e$INTRON.xlink)

d <- a[c("gid", "KD.TR", "Ctrl.TR", "KDdivCtrl.TR","INTRON.xlink","KD.INTRON","KD.aroundINTRON","aroundINTRON.length","INTRON.length")]
write.table(d,file="intron4TR.out",quote=F,sep="\t",row.names=F);

a[a$group=="CLIPPED",]

head(a,2)

xtabs(~  a$group)
for (g in unique(a$group)) {
	fivenum(a[a$group==g,]$KDdivCtrl.TR)
}
#ecdfplot(~ (KD.TR), group=group, data=a, 
ecdfplot(~ (Ctrl.TR), group=group, data=a, 
ecdfplot(~ KDdivCtrl.TR, group=group, data=a, 
	scales=list(x=list(log=T),y=list(log=F)),
	types=c("l","g"),
#	subset=length>1000,
	subset=aroundINTRON.length>300&INTRON.length>10000,
	auto.key=list(space="top"))

quit();
bwplot((KD.TR / Ctrl.TR) ~ group, data=a,
	panel = function(x, y, ...) {
	panel.grid(h = -1, v = -1)
	panel.bwplot(x, y, ...)
})


xyplot(KD.TR ~ Ctrl.TR, data=a,
	group=group,type=c("p","smooth","g"),pch=".",
#	subset=length>10000,
	auto.key=list(columns=2),
	scales=list(x=list(log=T),y=list(log=T)))

xyplot(displacementMedian ~ length, data=a,
	scales=list(x=list(log=T)),
	type=c("p","smooth"))


bwplot(displacementMedian + displacementSd, data=a)
fivenum(a$displacementMedian)
cor(a$displacementMedian,a$length, method="spearman")

