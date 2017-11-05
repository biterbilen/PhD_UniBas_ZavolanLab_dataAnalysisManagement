library(ggplot2)

df <- data.frame( 
	trt = factor(c(1, 1, 2, 2)), 
	resp = c(1, 5, 3, 4), 
	group = factor(c(1, 2, 1, 2)), 
	se = c(0.1, 0.3, 0.3, 0.2) 
) 
df2 <- df[c(1,3),] 

# Define the top and bottom of the errorbars 
limits <- aes(ymax = resp + se, ymin=resp - se) 

p <- ggplot(df, aes(fill=group, y=resp, x=trt)) 
p + geom_bar(position="dodge", stat="identity") 

dodge <- position_dodge(width=0.9) 
p + geom_bar(position=dodge) + geom_errorbar(limits, position=dodge, width=0.25) 

p <- ggplot(df, aes(colour=group, y=resp, x=trt)) 
p + geom_point() + geom_errorbar(limits, width=0.2) 

p + geom_crossbar(limits, width=0.2) 

p + geom_line(aes(group=group)) + geom_errorbar(limits, width=0.2) 

