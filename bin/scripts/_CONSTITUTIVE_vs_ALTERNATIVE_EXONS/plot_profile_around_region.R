cargs <- commandArgs();

library("latticeExtra") 

#cargs = c(1,1,1,"../Reproducibility_v06/_020_nhood150_topEvent10000000_Profile_.profile", "hnRNPC_SILENCED_CE");

a <- read.table(cargs[4], header=T, check.names=F)

eventName <- as.character(names(xtabs( ~eventName, data=a)))
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
	tit = paste("Event Profiles around ", cargs[5], sep="")
	tit
	ofilename <- paste(cargs[4], "_", cargs[5], "_nhood", snhood, ".pdf", sep="")
	plot.new(); pdf(ofilename, width=8, height=8)
	#xyplot(log2(eventCount) ~ distance | eventName, data=a,  #bad since the total number of events are different in each profile
	xyplot((100*eventCount/totalNumberOfEvents) ~ distance | eventName, data=a, #bad if u include x=0
	#xyplot(1/(-log(eventCount/totalNumberOfEvents)) ~ distance | eventName, data=a, 
	#xyplot(log(100 * 1000000 * eventCount/totalNumberOfEvents) ~ distance | eventName, data=a, 
		type=c("h","g"), 
		subset=region == cargs[5] &  distance > -(snhood+1) & distance < (snhood+1) & distance != 0 ,
		groups=eventName, auto.key=list(title=tit, columns=4), 
		layout=c(2,ceiling(length(eventName)/2))) 
	dev.off();
#}


