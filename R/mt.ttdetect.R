

#setwd("G:\\MPI\\__trajtypes_paper_2015\\RawData")
#data <- readRDS("koop_processed.RDS")
#library(ggplot2)

## input:
# data matrix with MAD column labeled "MAD"
# number of clusters: nclust
## output:
# data matrix with one additional column = cluster-membership


mt.ttdetect <- function(data, nclust) {
  
fit <- kmeans(data[,c("MAD")], nclust)
data$clusters <- fit$cluster

return(data)

} # end of function


#table(data$clusters)


#visualize
#setwd("G:\\MPI\\__trajtypes_paper_2015\\analysis")

#pdf(file="kmeans2.pdf", height=8, width=11)
#p <- ggplot(data, aes(xflip,y))
#p + geom_point(alpha=.5, color="#FF6666") + theme_bw()
#p <- ggplot(data, aes(xflip,y))
#p + geom_point(aes(colour = factor(cluster)), alpha=.5) + theme_bw() + theme(legend.position="none")
#dev.off()

#boxplot(MAD~cluster, data)





