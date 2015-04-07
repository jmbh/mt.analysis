## input:
# 1. mouse tracking data
# 2. indicating of vector we should use for aggregation

## output:
# 1. (1) aggregated with respect to specified vector

#library(plyr)

mt.aggregate <- function(
  data, #mousetracking data
  i.aggr, #variable using for aggregation
  i.xyt #x,y  and time values
) {
  
out2 <- ddply(data, i.aggr, function(z) {
  
  out1 <- ddply(z, c("t"), function(x) {
   
   mx <- mean(as.matrix(x[i.xyt[1]]))
   my <- mean(as.matrix(x[i.xyt[2]]))
   
  return(cbind(mx, my))
    })
  
  return(out1)

  })

return(out2)

} #end of function
  

## TESTDATA
#setwd("G:\\MPI\\__trajtypes_paper_2015\\RawData")
#data <- readRDS("koop_processed.RDS")
#i.aggr <- "hrisk"
#i.xyt <- c("xflip", "y", "t")
#agg <- mt.aggregate(data, i.aggr, i.xyt)

