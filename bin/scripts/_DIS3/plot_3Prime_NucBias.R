library(lattice)

cargs <- commandArgs();
#cargs <- c(1,1,1,"woG0");
#cargs <- c(1,1,1,"woG1.fpA.enrch");

a <- read.table(cargs[4], header=T);

ltheme <- canonical.theme(color = F)     ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme)      ## set as default

#trellis.par.set(canonical.theme(color = FALSE))
#trellis.par.set(col.whitebg()) 

plot.new(); pdf(paste(cargs[4],'.pdf',sep=''),h=8,w=8)
#barchart(A/posNucSum+C/posNucSum+G/posNucSum+T/posNucSum~ pos | factor(lib), data=a, 
barchart(log2(Aenrch) + log2(Cenrch) + log2(Genrch) + log2(Tenrch) ~ pos | factor(lib), data=a, 
#barchart(Afreq + Cfreq + Gfreq + Tfreq ~ pos | factor(lib), data=a, 
	origin=0, layout=c(1,7),auto.key=list(space="right"), main=cargs[4]);
dev.off();

