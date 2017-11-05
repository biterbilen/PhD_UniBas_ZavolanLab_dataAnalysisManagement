cargs <- commandArgs();
library(lattice)

#######
#cargs = c(1,1,1,"Reproducibility_v05T2C_Profile_.stat")
a <- read.table(cargs[4], header=T, check.names=F)

nhood <- a$nhood[1];
totalEventCount <- a$totalNumberOfEvents[1];
totalRegionCount <- a$totalNumberOfRegions[1];
totalHitRegionCount <- a$totalNumberOfHitRegions[1];

tit = paste("Adjusted Regional Preference of Event Profiles", sep="")
plot.new(); pdf(paste(cargs[4],'.pdf',sep=''),w=8,h=8)
barchart(~totalNumberofHitRegions/totalNumberOfRegions | eventName, data=a,
	type=c('g'), groups=region, stack=T,
	auto.key=list(title=tit, columns=2))
dev.off()

quit();
#-------------------------------------------
cargs = c(1,1,1,"Reproducibility_v05ASeq_Profile.stat")

a <- read.table(cargs[4], header=T, check.names=F)

nhood <- a$nhood[1];
totalEventCount <- a$totalNumberOfEvents[1];
totalRegionCount <- a$totalNumberOfRegions[1];
totalHitRegionCount <- a$totalNumberOfHitRegions[1];

tit = paste("Regional Preference of Top ", totalEventCount," ASEQ Profiles", sep="")

plot.new(); pdf(paste(cargs[4],'.pdf',sep=''),w=8,h=8)
barchart(~totalNumberofHitRegions | eventName, data=a,
	type=c('g'), groups=region, stack=T,
	auto.key=list(title=tit, columns=3))
dev.off()

