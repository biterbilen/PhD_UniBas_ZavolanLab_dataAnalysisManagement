set.seed(100)
Data <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
	matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
colnames(Data) <- c("x", "y")

# noise <- runif(100,0,.05)
line.width <- rep(.004, dim(Data)[1])
Y <- NULL
X <- NULL
k.range <- 2:10

plot(0,0, col = "white", xlim = c(1,10), ylim = c(-.5,1.6),
	xlab = "Number of clusters", ylab = "Clusters means", main = "(Basic) Clustergram")
axis(side =1, at = k.range)
abline(v = k.range, col = "grey")

centers.points <- list()

for(k in k.range)
{
	cl <- kmeans(Data, k)

	clusters.vec <- cl$cluster
	the.centers <- apply(cl$centers,1, mean)

	noise <- unlist(tapply(line.width, clusters.vec, cumsum))[order(seq_along(clusters.vec)[order(clusters.vec)])]  
	noise <- noise - mean(range(noise))
	y <- the.centers[clusters.vec] + noise
	Y <- cbind(Y, y)
	x <- rep(k, length(y))
	X <- cbind(X, x)

	centers.points[[k]] <- data.frame(y = the.centers , x = rep(k , k)) 
	#   points(the.centers ~ rep(k , k), pch = 19, co
}	

require(colorspace)
COL <- rainbow <- hcl(100)
matlines(t(X), t(Y), pch = 19, col = COL, lty = 1, lwd = 1.5)

# add points
lapply(centers.points, function(xx) {with(xx,points(y~x, pch = 19, col = "red", cex = 1.3))})

