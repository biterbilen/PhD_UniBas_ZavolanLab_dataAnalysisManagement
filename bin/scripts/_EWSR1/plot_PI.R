cargs <- commandArgs()

#cargs <- c(1,1,1,"Expression/RNAseq.raw.qn.w_geneSymbol.FC.PI.TF");

a <- read.table(cargs[4], header=T)
head(a,2)

library(latticeExtra)
p1 <- ecdfplot(~mRNAseq_EWSR1_siEWSR1_1_over_mRNAseq_EWSR1_siCTRL_1, groups=PIstat, data=a,
	main="RNAseq FC for PIstat",
	scales = list(x = list(rot = 90)),
	auto.key=list(space = "top", title = "EWSR1 interaction") )

p2 <- ecdfplot(~mRNAseq_EWSR1_siEWSR1_2_over_mRNAseq_EWSR1_siCTRL_2, groups=PIstat, data=a,
	scales = list(x = list(rot = 90)),
	auto.key=list(space = "top", title = "EWSR1 interaction") )

plot.new(); pdf(paste(cargs[4],"_foldchange.pdf",sep=""),h=8,w=8);
plot(p1, split = c(1, 1, 1, 2)) 
plot(p2, split = c(1, 2, 1, 2), newpage = FALSE) 
dev.off()

alt <- "two.sided";
alt <- "greater"
s <- "batch1"
xi <- a$PIstat=="PI";
x <- a[xi,"mRNAseq_EWSR1_siEWSR1_1_over_mRNAseq_EWSR1_siCTRL_1"];
yi <- a$PIstat=="NOTPI";
y <- a[yi,"mRNAseq_EWSR1_siEWSR1_1_over_mRNAseq_EWSR1_siCTRL_1"];
ks <- ks.test(x,y,alternative=alt)
write.table(data.frame(list(method="Two sample KS test for EWSR1 interacting and not interacting RNAseq expression FC", alternative=alt, sampleName=s, pvalue=ks$p.value, statistic=ks$statistic)),sep="\t",row.names=F, quote=F);

s <- "batch2"
x <- a[xi,"mRNAseq_EWSR1_siEWSR1_2_over_mRNAseq_EWSR1_siCTRL_2"];
y <- a[yi,"mRNAseq_EWSR1_siEWSR1_2_over_mRNAseq_EWSR1_siCTRL_2"];
ks <- ks.test(x,y,alternative=alt)
write.table(data.frame(list(method="Two sample KS test for EWSR1 interacting and not interacting RNAseq expression FC", alternative=alt, sampleName=s, pvalue=ks$p.value, statistic=ks$statistic)),sep="\t",row.names=F, quote=F);


#----------------------------------
p1 <- ecdfplot(~mRNAseq_EWSR1_siEWSR1_1_over_mRNAseq_EWSR1_siCTRL_1, groups=PIstat, data=a,
	subset=!is.na(TF),
	main="TF RNAseq FC for PIstat",
	scales = list(x = list(rot = 90)),
	auto.key=list(space = "top", title = "EWSR1 interaction") )

p2 <- ecdfplot(~mRNAseq_EWSR1_siEWSR1_2_over_mRNAseq_EWSR1_siCTRL_2, groups=PIstat, data=a,
	subset=!is.na(TF),
	scales = list(x = list(rot = 90)),
	auto.key=list(space = "top", title = "EWSR1 interaction") )

plot.new(); pdf(paste(cargs[4],"_TF_foldchange.pdf",sep=""),h=8,w=8);
plot(p1, split = c(1, 1, 1, 2)) 
plot(p2, split = c(1, 2, 1, 2), newpage = FALSE) 
dev.off()

a
plot.new(); pdf(paste(cargs[4],"_activity.pdf",sep=""),h=8,w=8);
ecdfplot(~zvalue, groups=PIstat, data=a,
	subset=!is.na(TF),
	main="TF activiy for PIstat",
	scales = list(x = list(rot = 90)),
	auto.key=list(space = "top", title = "EWSR1 interaction") )
dev.off()


alt <- "two.sided"
alt <- "greater"
s <- "batch1"
xi <- a$PIstat=="PI"&!is.na(a$TF);
x <- a[xi,"mRNAseq_EWSR1_siEWSR1_1_over_mRNAseq_EWSR1_siCTRL_1"];
yi <- a$PIstat=="NOTPI"&!is.na(a$TF);
y <- a[yi,"mRNAseq_EWSR1_siEWSR1_1_over_mRNAseq_EWSR1_siCTRL_1"];
ks <- ks.test(x,y,alternative=alt)
write.table(data.frame(list(method="Two sample KS test for EWSR1 interacting and not interacting TFs RNAseq expression FC", alternative=alt, sampleName=s, pvalue=ks$p.value, statistic=ks$statistic)),sep="\t",row.names=F, quote=F);

s <- "batch2"
x <- a[xi,"mRNAseq_EWSR1_siEWSR1_2_over_mRNAseq_EWSR1_siCTRL_2"];
y <- a[yi,"mRNAseq_EWSR1_siEWSR1_2_over_mRNAseq_EWSR1_siCTRL_2"];
ks <- ks.test(x,y,alternative=alt)
write.table(data.frame(list(method="Two sample KS test for EWSR1 interacting and not interacting TFs RNAseq expression FC", alternative=alt, sampleName=s, pvalue=ks$p.value, statistic=ks$statistic)),sep="\t",row.names=F, quote=F);

#alt <- "less"
s <- "batch12"
x <- a[xi,"zvalue"];
y <- a[yi,"zvalue"];
ks <- ks.test(x,y,alternative=alt)
write.table(data.frame(list(method="Two sample KS test for EWSR1 interacting and not interacting TFs MARA activity", alternative=alt, sampleName=s, pvalue=ks$p.value, statistic=ks$statistic)),sep="\t",row.names=F, quote=F);

