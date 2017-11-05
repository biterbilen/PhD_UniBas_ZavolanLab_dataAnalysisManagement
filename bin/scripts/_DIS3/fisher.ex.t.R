#PAPD5 parameters
#a <- matrix(c(2,29,30,71), nr=2); #before
#fisher.test(a, alternative='greater')
cargs <- commandArgs();
#cargs <- c(1,1,1,"../Normalization_PAPD5_T2C/_04_14.csv.1.stat");
a <- read.table(cargs[4], header=F);

v <- "Total";

dd <- c();
for (v2 in c("Total", "Repeat")) {
	totalEnriched <- a[a$V2=="Enriched" & a$V1==v2,3];
	totalDepleted <- a[a$V2=="Depleted" & a$V1==v2,3];
	#v<-"Repeat_and_rRNA"
	for (v in levels(a$V1)) {
		if (v == "Total") { next; }
		#skip Repeat over Repeat comparison
		if ( v == v2 && v == "Repeat") { next; } 
		#skip over Repeat comparison if the category is not a repeat
		if (v2 == "Repeat" && length(grep(v2,v,value=T)) < length(v2)) { next; }
		#write.table(c(v, v2));
		ep <- a[a$V2=="Enriched" & a$V1==v,3];
		dp <- a[a$V2=="Depleted" & a$V1==v,3];
		d <- matrix(c(ep,dp,totalEnriched-ep,totalDepleted-dp), nr=2); 
		alt <- "greater"
		if (dp > ep ) {
			alt <- "less";
		}
		ts <- fisher.test(d, alternative=alt);
		#print(ts);
		dd <- cbind(dd, c(v, paste(alt,"_o_",v2, sep=""), as.vector(d), ts$p.value));
		#print(dd);
	}
}
write.table(t(dd), file=paste(cargs[4], '.pval', sep=""), quote=F, col.names=F, row.names=F);

