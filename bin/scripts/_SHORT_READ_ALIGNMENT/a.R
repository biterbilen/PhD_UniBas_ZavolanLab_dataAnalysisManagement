library(lattice)

cargs <- c(1,1,1,"transcripts.EXP");
a <- read.table(cargs[4], header=T)
hall <- histogram(~log2(score),data=a,nint=50, type="count", subset=score>100,main="all")
hpart <- histogram(~log2(score),data=a,nint=50, type="count", subset=score>100 && a.exp>0 && b.exp>0, main="score>100"))

plot.new(); pdf(paste(cargs[4],".score_histogram.pdf",sep=""),h=8,w=8);
plot(hall, split = c(1, 1, 1, 2));
plot(hpart,split = c(1, 2, 1, 2), newpage = FALSE)
dev.off()

