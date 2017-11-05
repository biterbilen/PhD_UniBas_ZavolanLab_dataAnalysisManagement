library(latticeExtra)
library(hexbin)

ltheme <- canonical.theme(color = F)     ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme)      ## set as default

cargs <- c(1,1,1,"../Normalization_PAPD5/_01_01.filtered");


cargs <- c(1,1,1,"../Normalization_PAPD5/tableall.norm_and_raw.01zscore.woNA.filtered");
a <- read.table(cargs[4], header=F)
d01 <- data.frame(comp="01", zscore=a$V6);

cargs <- c(1,1,1,"../Normalization_PAPD5/tableall.norm_and_raw.04zscore.woNA.filtered");
a <- read.table(cargs[4], header=F)
d04 <- data.frame(comp="04", zscore=a$V6);

cargs <- c(1,1,1,"../Normalization_PAPD5/tableall.norm_and_raw.14zscore.woNA.filtered");
a <- read.table(cargs[4], header=F)
d14 <- data.frame(comp="14", zscore=a$V6);

my.data <- rbind(d01, d04, d14);
my.filt <- my.data[my.data$zscore>4,];
my.filt <- my.data;
#ecdfplot(~ zscore | equal.count(my.filt$zscore, 1, overlap=0), data = my.filt, groups = comp, auto.key = list(columns = 2), xlab = "Average zScore")
ecdfplot(~ zscore | cut(my.filt$zscore, breaks = c(-10,-6,-3,3,6,10)), data = my.filt, groups = comp, auto.key = list(columns = 2), xlab = "Average zScore")
#ecdfplot(~ zscore | cut(my.filt$zscore, 4), data = my.filt, groups = comp, auto.key = list(columns = 2), xlab = "Average zScore")


cargs <- c(1,1,1,"../Normalization_PAPD5/_04_14.filtered");
a <- read.table(cargs[4], header=F, col.names=c("id", "x1.e", "y.e", "x1y.z", "x2.e", "y.e", "x2y.z"));

xyplot(log(x1.e) * log(x2.e) ~ log(y.e), data= a, col=1, type=c("g","p"), panel=panel.hexbinplot) 

k <- xyplot(log(x1.e) ~ log(y.e), data= a, col=1, type=c("g","p"), panel=panel.hexbinplot) +
#k1 <- xyplot(log(x2.e) ~ log(y.e), data= a, col=1, type=c("g","p"), panel=panel.hexbinplot) +
layer(panel.abline(a=0,b=1, col=1)) +
layer(panel.key(c(sprintf('R=%.3f N=%d', cor(x,y), length(x))),lwd=2,cex=2)) +
as.layer(xyplot(log(x1.e) ~ log(y.e), data= a, subset=x1y.z > 3, pch=1)) +
as.layer(xyplot(log(x1.e) ~ log(y.e), data= a, subset=x1y.z < -3, pch=1)) +
as.layer(xyplot(log(x1.e) ~ log(y.e), data= a, subset=x1y.z > 4, pch=2)) +
as.layer(xyplot(log(x1.e) ~ log(y.e), data= a, subset=x1y.z < -4, pch=2));

plot.new(); pdf(paste(cargs[4], '.pdf', sep=""), h=8,w=8);
update(c(k, k1), layout=c(1,2));
dev.off();

