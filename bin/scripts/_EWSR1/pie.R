library(lattice)
source(file.path("/import/bc2/home/zavolan/bilebi00/_ARE/scripts/","theme.R"))

cargs <- commandArgs(trailingOnly=T);
#cargs <- c("../Summary/all.mRNA");
#cargs <- c("annot.stats")
a <- read.table(cargs[1], header=T, sep="\t");
#sum reads
b <- aggregate(value ~ lib, a, sum)
names(b) <- gsub("value", "sum", names(b))
#get freq
b <- merge(a,b, by=c("lib"))
b$frequency <- b$value / b$sum
#mark small freqs to remove
b$other <- b$frequency < 0.01
#summ small freqs
d <- aggregate(value ~ lib, b[b$other,], sum)
d$type <- "other"
#merge
e <- rbind(b[b$other==F,names(d)],d)
head(e)
#prep input for pie
b <- xtabs(value ~ type + lib, e, drop.unused.levels = T)

panel.piechart <-
function(x, y, labels = as.character(y),
	edges = 200, radius = 0.8, clockwise = FALSE,
	init.angle = if(clockwise) 90 else 0,
	density = NULL, angle = 45, 
	col = superpose.polygon$col,
	border = superpose.polygon$border,
	lty = superpose.polygon$lty, ...)
{
	stopifnot(require("gridBase"))
	superpose.polygon <- trellis.par.get("superpose.polygon")
	opar <- par(no.readonly = TRUE)
	on.exit(par(opar))
	if (panel.number() > 1) par(new = TRUE)
	par(fig = gridFIG(), omi = c(0, 0, 0, 0), mai = c(0, 0, 0, 0))
	pie(as.numeric(x), labels = cbind(labels, x), edges = edges, radius = radius,
#	pie(as.numeric(x), labels = labels, edges = edges, radius = radius,
		clockwise = clockwise, init.angle = init.angle, angle = angle,
		density = density, col = col, border  = border, lty = lty)
}

piechart <- function(x, data = NULL, panel = "panel.piechart", ...)
{
	ocall <- sys.call(sys.parent())
	ocall[[1]] <- quote(piechart)
	ccall <- match.call()
	ccall$data <- data
	ccall$panel <- panel
	ccall$default.scales <- list(draw = FALSE)
	ccall[[1]] <- quote(lattice::barchart)
	ans <- eval.parent(ccall)
	ans$call <- ocall
	ans
}

#order
nms <- names(sort(apply(b,1,sum), decreasing=T))
b <- b[nms,]

plot.new(); pdf(paste(cargs[1], '.pdf', sep=''), h=8, w=8);
#data(NHANES) 
#piechart(VADeaths, groups = FALSE, xlab = "") 
ltheme <- canonical.theme(color = FALSE)      ## in-built B&W theme
ltheme$strip.background$col <- "transparent" ## change strip bg
lattice.options(default.theme = ltheme)      ## set as default
#lattice.options(default.theme = article.theme)      ## set as default
#trellis.device(color = FALSE)
piechart(b, groups=FALSE, layout=c(1,2))
dev.off();

