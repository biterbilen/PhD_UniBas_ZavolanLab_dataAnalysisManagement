cargs <- commandArgs();
library(lattice)

#input is an expression matrix where the rows are the uniq test id, values, gene symbol respectively
#assuming multiplicative counts [each entry in a row + pseoudocount of 1] is log2 transformed and divided by the row maximum 
#which is similar to codon adaptation index where the denominator is the log2(row maximum) 
#i.e. each row is scaled to 0-1
#
#cargs <- c(1,1,1,"../Stepanka/Parallel_plots/dima_clip_targets_in_Hela.clipzQNgene");
a <- read.table(cargs[4], header=T);
rs <- dim(a)[1]
cs <- dim(a)[2] - 1;
d <- log2(a[,2:cs] + 1);
m <- apply(d,1,max)
M <- matrix(rep(m, cs),c(rs,cs));
D <- data.frame(a$id,a$gid, d/M)
#L <- reshape(D, timevar="id", idvar="gid", varying=list(names(d)), direction="long");
L <- reshape(D, timevar="id", idvar="gid", varying=list(names(d)), direction="long",v.names="expr");

ltheme <- canonical.theme(color = FALSE)      ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme, lwd=4)      ## set as default

plot.new(); pdf(paste(cargs[4],'.01scaled.parallel.pdf', sep=''),h=8,w=8);
xyplot(expr ~ id | a.gid, data=L, type=c("g","l"),
	layout=c(3,1),
	par.settings=list(axis.text=list(font=2,cex=1.0),
		par.strip.text=list(font=2,cex=1.0),
		par.ylab.text=list(font=2,cex=1.0),
		par.xlab.text=list(font=2,cex=1.0),
		par.main.text=list(font=2,cex=1.0)),
	lwd=4,
	scales = list(x = list(rot=90,at=1:cs, labels=names(d))),
	xlab="", main=cargs[4], ylab="Scaled Log2 Expression")
dev.off();

#previous tries where the mx and min resolution was a problem
#plot.new(); pdf(paste(cargs[4],'.parallel.pdf', sep=''),h=8,w=8);
#parallel(~(a[,2:cs])|gid, data=a, xlab="Log Expression Index",common.scale=F, main=cargs[4],
##	scales = list(x = list(rot=90, at = (prng-min(rng))/diff(rng), labels = prng)));
#	scales = list(x = list(rot=90)));
#dev.off();

#plot.new(); pdf(paste(cargs[4],'.parallel.pdf', sep=''),h=8,w=8);
#parallel(~ asinh(a[c(2:12)]), data = a, alpha = 0.1, lty = 1) 
#dev.off()

