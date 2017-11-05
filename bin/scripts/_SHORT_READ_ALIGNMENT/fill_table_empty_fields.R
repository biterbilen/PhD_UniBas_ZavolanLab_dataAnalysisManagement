cargs <- commandArgs();

cargs <- c(1,1,1,"../data/Splice_arrays/Harald_31_03_2011_mouse_cerebellum/splicing-change-Harald_31_03_2011_mouse_cerebellum-High_Low_with_duplicates.txt");

a <- read.table(cargs[4], header=T, sep="\t",  na.strings="NA")

write.table(a, file=paste(cargs[4],".iinp", sep=""), quote=F, sep="\t", row.names=F) 

