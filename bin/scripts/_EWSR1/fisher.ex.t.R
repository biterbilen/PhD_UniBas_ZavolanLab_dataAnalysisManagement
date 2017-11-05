cargs <- commandArgs(trailingOnly=T);
#library(qvalue)

#cargs <- c(29,799,229,21336);

clippedChimer <- as.numeric(cargs[1]);
clipped <- as.numeric(cargs[2]);
chimer <- as.numeric(cargs[3]);
total <- as.numeric(cargs[4]);
cont <- 0;
if (length(cargs) > 4) {
	cont <- as.numeric(cargs[5]); #the values are for contingency table
}

paste(clipped)

d <- matrix(c(clippedChimer,clipped,chimer, total), nr=2); 
if (cont == 0) {
	d <- matrix(c(clippedChimer,clipped-clippedChimer,chimer-clippedChimer, total-(clipped+chimer-clippedChimer)), nr=2); 
}
paste(d)

alt <- "greater"

tt <- fisher.test(d, alternative=alt);

if (cont == 0) {
	write.table(format(digits=2, data.frame(list(method="Fisher.test", alternative=alt, set12count=clippedChimer, set1count=clipped, set2count=chimer, totalcount=total, pvalue=tt$p.value, oddsratio=tt$estimate))),sep="\t",row.names=F, quote=F);
} else {
	write.table(format(digits=2, data.frame(list(method="Fisher.test", alternative=alt, X1Y1=clippedChimer, X1Y2=clipped, X2Y1=chimer, X2Y2=total, pvalue=tt$p.value, oddsratio=tt$estimate))),sep="\t",row.names=F, quote=F);
}


