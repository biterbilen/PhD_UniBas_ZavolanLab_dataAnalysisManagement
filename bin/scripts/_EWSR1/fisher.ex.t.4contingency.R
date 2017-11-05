cargs <- commandArgs(trailingOnly=T);

#cargs <- c(29,799,229,21336);

clippedChimer <- as.numeric(cargs[1]);
clipped <- as.numeric(cargs[2]);
chimer <- as.numeric(cargs[3]);
total <- as.numeric(cargs[4]);

paste(clipped)

d <- matrix(c(clippedChimer,clipped-clippedChimer,chimer-clippedChimer, total-(clipped+chimer-clippedChimer)), nr=2); 
paste(d)

alt <- "greater"

tt <- fisher.test(d, alternative=alt);

write.table(format(digits=2, data.frame(list(method="Fisher.test", alternative=alt, set12count=clippedChimer, set1count=clipped, set2count=chimer, totalcount=total, pvalue=tt$p.value, oddsratio=tt$estimate))),sep="\t",row.names=F, quote=F);


