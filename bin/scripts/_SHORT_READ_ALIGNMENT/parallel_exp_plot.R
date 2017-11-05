cargs <- commandArgs();
library(lattice)

#cargs <- c(1,1,1,"~bilebi00/_SHORT_READ_ALIGNMENT/Stepanka/Parallel_plots/clipzQNgene.dima_clip_targets_in_Hela");
#cargs <- c(1,1,1,"~bilebi00/_SHORT_READ_ALIGNMENT/Stepanka/Parallel_plots/CuffdiffTrustedgene.dima_clip_targets_in_Hela");
cargs <- c(1,1,1,"~bilebi00/_SHORT_READ_ALIGNMENT/Stepanka/Parallel_plots/dima_clip_targets_in_Hela.CuffdiffTrustedisoform");
cargs <- c(1,1,1,"~bilebi00/_SHORT_READ_ALIGNMENT/Stepanka/Parallel_plots/DIS3s.CuffdiffTrustedgene");
cargs <- c(1,1,1,"~bilebi00/_SHORT_READ_ALIGNMENT/Stepanka/Parallel_plots/DIS3s.clipzQNgene");
a <- read.table(cargs[4], header=T);
cs <- dim(a)[2] - 1;

b <- a[,2:cs];

rng <- range(unlist(lapply(lapply(b, as.numeric), range)))
prng <- pretty(rng)


plot.new(); pdf(paste(cargs[4],'.parallel.pdf', sep=''),h=8,w=8);
parallel(~(a[,2:cs])|gid, data=a, xlab="Scaled Expression",common.scale=F, main=cargs[4],
#	scales = list(x = list(rot=90, at = (prng-min(rng))/diff(rng), labels = prng)));
	scales = list(x = list(rot=90)));
dev.off();

#plot.new(); pdf(paste(cargs[4],'.parallel.pdf', sep=''),h=8,w=8);
#parallel(~ asinh(a[c(2:12)]), data = a, alpha = 0.1, lty = 1) 
#dev.off()

