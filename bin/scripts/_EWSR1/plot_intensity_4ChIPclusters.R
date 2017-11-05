cargs <- commandArgs(trailingOnly=T)

#cargs <- c("1.TfbsClustered", "CTCF", 1)

a <- read.table(cargs[1], header=F, sep="\t")
names(a) <- c("gene_id", "TF_name", "clipped_length", "clipped_intensity", "not_clipped_length", "not_clipped_intensity")

alt <- "greater"
tt <- t.test(a$clipped_intensity / a$clipped_length, a$not_clipped_intensity / a$not_clipped_length,alt=alt)
write.table(format(digits=2, data.frame(list(method=paste("Two sample t-test for average TF binding intensity (",cargs[2],")",sep=''), alternative=alt, sample1="clipped", sample2="notclipped", pvalue=tt$p.value, statistic=tt$statistic))),sep="\t",row.names=F, quote=F);

if (tt$p.value < 1e-1 | cargs[3] == 1) {
	library(latticeExtra)
	library(hexbin)

	b <- rbind(data.frame(type="clipped", averageIntensity=a$clipped_intensity / a$clipped_length), 
		data.frame(type="notclipped", averageIntensity=a$not_clipped_intensity / a$not_clipped_length));
	print(head(b))

	plot.new(); pdf(paste(cargs[2], '.pdf', sep=''), h=8, w=8) 
	print(densityplot( ~log2(averageIntensity), groups=type, data = b, plot.points=F, auto.key=list(corner=c(1,1)), xlab = paste("average intensity of ", cargs[2], sep='')));
	dev.off()
}
