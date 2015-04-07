

setwd("G:\\MPI\\__trajtypes_paper_2015\\RawData")
data <- readRDS("koop_processed.RDS")
library(ggplot2)

names(data)

#k-means
hist(data$MAD)
?kmeans
fit <- kmeans(data[,c("MAD")], 2)
fit

cluster <- fit$cluster
data$clusters <- cluster
table(data$clusters)

boxplot(MAD~cluster, data)


#visualize
setwd("G:\\MPI\\__trajtypes_paper_2015\\analysis")

#pdf(file="kmeans2.pdf", height=8, width=11)
p <- ggplot(data, aes(xflip,y))
p + geom_point(alpha=.5, color="#FF6666") + theme_bw()
p <- ggplot(data, aes(xflip,y))
p + geom_point(aes(colour = factor(cluster)), alpha=.5) + theme_bw() + theme(legend.position="none")
#dev.off()



#get lines instead....

dt <- data[1:101,]

p <- ggplot(data, aes(xflip,y))
p + geom_point(alpha=1, color="#FF6666") + theme_bw()






#
data1 <- subset(data, cluster==1)
data2 <- subset(data, cluster==2)
data3 <- subset(data, cluster==3)

p <- ggplot(data2, aes(xflip,y))
p + geom_point(alpha=.1)







