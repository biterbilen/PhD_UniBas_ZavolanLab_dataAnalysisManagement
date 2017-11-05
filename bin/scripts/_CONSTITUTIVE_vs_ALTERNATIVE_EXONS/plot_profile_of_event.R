cargs <- commandArgs();

library("latticeExtra") 

#cargs = c(1,1,1,"../Reproducibility_v06/_020_nhood150_topEvent10000000_Profile_.profile", "T1_hnRNPC");

a <- read.table(cargs[4], header=T, check.names=F)

#eventName <- as.character(names(xtabs( ~eventName, data=a)))
region <- as.character(names(xtabs( ~region, data=a)))
nhood <- a$nhood[1];
totalEventCount <- a$totalNumberOfEvents[1];
totalRegionCount <- a$totalNumberOfRegions[1];
totalHitRegionCount <- a$totalNumberOfHitRegions[1];

#Loop doesn't work maybe it can be given as a parameter to the program 
#for (f in c(1,4,16)) { 
#	f=4
#	snhood <- round(nhood / as.numeric(cargs[6]));
	snhood=nhood
	#tit = paste("Top ", totalEventCount," Event Profiles around ", cargs[5], sep="")
	tit = paste(cargs[5], " Profiles around Regions", sep="")
	ofilename <- paste(cargs[4], "_", cargs[5], "_nhood", snhood, ".pdf", sep="")
	ofilename
	plot.new(); pdf(ofilename, width=8, height=8)

	xyplot((100*eventCount/totalNumberOfEvents) ~ distance | region, data=a, #bad if u include x=0
		type=c("h","g"), 
		subset= eventName == cargs[5] & distance > -(snhood+1) & distance < (snhood+1) & distance != 0,
		groups=region, auto.key=list(title=tit, columns=2),
		layout=c(2,ceiling(length(region)/2))) 
	dev.off();
#}


