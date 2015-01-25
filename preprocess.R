####################################################
### Preprocess mouse-tracking data for analysis
### jonashaslbeck@gmail.com
### January 2015
####################################################

#depedencies
library(plyr)

## input:
# 1) data matrix containing the columns: trialid (unique id for each trial/trajectory),x,y,time,chosen box (left/right)
# 2) coordinates of the starting position and the two boxes of the mouse

box.cor <- list("start"=c(960,230), "left"=c(130,905), "right"=c(1830,905)) #coordinates

## output:
# 1) time-normalized x and y values lined up to (0,0); x also fliped on one side
# 2) on time step level: velocity(difference), distance to direct line
# 3) on trial level: RT, mean velocity, MAD, AAD, directionalchange(DC), total distance


#TEST DATA

setwd("G:/MPI/_trajtypes_paper_2015")
load("testdata_experience.RData")
setwd("G:/MPI/mt.analysis")

data <- d_test <- test[,c("trial_count", "x", "y", "time", "box_chosen")]
rownames(data) <- NULL

box.cor <- list("start"=c(960,230), "left"=c(130,905), "right"=c(1830,905)) #coordinates


mt.preprocess <- function(
  data, #see above
  box.cor #see above
  ) {
  

  
  #### step 1 - sanity checks 

  if(length(box.cor) != 3) {
    stop("Please specify the positions of the boxes.")
    }
  if(ncol(data)!=5) {
    stop("Please provide the data matrix in the right format.")
    }
  if(sum(is.na(data)>0)) {
    stop("No missing values allowed.")
    }

data <- as.data.frame(data)
colnames(data) <- c("id", "x", "y", "t", "b")
  

  #### step 2 -  set starting point to zero
  
  f_settozero <- function(z) 
  {
    x <- z[1] - box.cor$start[1]
    y <- z[2] - box.cor$start[2]
    return(rbind(x,y))
  }
  
  new.cord <- t(apply(data[,2:3], 1, f_settozero))
  data$x <- new.cord[,1]
  data$y <- new.cord[,2]
  

  #### step 3 - calc reaction time
  
  f_calc_rts <- function(x) 
  {
    le <- length(x$t)
    rt <- x$t[le] - x$t[1]
    return(cbind(rep(rt, le)))
  }
  x <- data

  data$rt <- ddply(data, c("id"), f_calc_rts)[,2]
  

  #### step 4 - time normalise ####
   

  f_timenorm <- function(x) 
  {
    v_x <- x[,2]
    v_y <- x[,3]
    v_time <- x[,4]
    
    #normalise time vector [0,1]
    v_time_norm <- (v_time - v_time[1]) / (v_time[length(v_time)]-v_time[1])
    v_time_norm <- v_time_norm*100
    
    # interpolation of x&y to 101 equally spaced time slices
    lin.x <- approx(v_time_norm, v_x, xout = 0:100, method = "linear") # interpolate x coordinates
    lin.y <- approx(v_time_norm, v_y, xout = 0:100, method = "linear") # interpolate y coordinates
    
    rawdata_restoftable1 <- x[rep(1,101),1] 
    rawdata_restoftable2 <- x[rep(1,101),5:6]

    
    data_export <- cbind(rawdata_restoftable1, lin.x$y, lin.y$y, lin.x$x, rawdata_restoftable2) 
        return(data_export)
  }

  
  data <- ddply(data, c("id"), f_timenorm)[,-1]
  colnames(data) <- c("id", "x", "y", "t", "b", "rt")


  #### step 5 - flip right-trajs to the left

  #flip all trajectories to the left side (arbitrary choice)

  f_flip <- function(z)
  {
    x <- z[2]
    if(z[1] == 1) 
    {
      x <- x * (-1)
    } 
    return(x)
  }
  
  xflip <- apply(cbind(data[,5], data[,2]),1,f_flip)
  data$xflip <- xflip



  #### step 6 - M: distance to direct line ####
  
  f_calculatedistance<- function(z)
  {
    
    slopedl <- (box.cor$left[2]-box.cor$start[2])/(box.cor$left[1]-box.cor$start[1]) #slope of direct line
    
    #find line with slope = -1/slopedl, which passes through the empirical point (interesting thing is only the intercept with y-axis, because slopes are known)
    c_x<-z[7]
    c_y<-z[3]
    slopedl_perp <- -1/slopedl
    c <- - slopedl_perp * c_x + c_y
    
    #calc x of intersection point; x = (-intercept1 + intercept2) / slope1 - slope2 (simple solving of equation of both linear equations)
    x_intersect <- (-0 + c) / (slopedl-slopedl_perp)
    y_intersect <- x_intersect * slopedl
    #sanity <- x_intersect * slopedl_perp + c # should be the same
    
    #euclidian distance between both points
    c_distance <- sqrt((x_intersect - c_x)^2 + (y_intersect - c_y)^2)
    return(c_distance)
  }
  
  data$dist <- apply(data, 1, f_calculatedistance) # raw-data
  

  #### step 7 - M: AAD/MAD #### 
  
  f_aadmad <- function(x) {
    AAD <- rep(mean(x$dist),101) # aad
    MAD <- rep(max(x$dist),101) # mad
    out <- cbind(AAD, MAD)
    return(out)
  }
  
  v_aadmad <- ddply(data, c("id"), f_aadmad)
  data$AAD <- v_aadmad$AAD
  data$MAD <- v_aadmad$MAD
  
  
  #### step 8 - M: Velocity ####
  
  dx <- c(0,data$xflip[-length(data$xflip)]) #dummy variables, vector shifted by one
  dy <- c(0,data$y[-length(data$y)])
  velo <- sqrt((data$xflip-dx)^2+(data$y-dy)^2) #euclidean distance between (x,y) at t-1 and t
  
  
  data$velo <- velo
  data$velo[data$t==0] <- NA #there cant be velocity in the first timestamp (no reference before)
  
  
  #### step 9 - M: Total distance / mean Velocity ####
    sumdist <- ddply(data, c("id"), function(x) {
    totdist <- sum(x$velo, na.rm=TRUE)
    meanvelo <- mean(x$velo, na.rm=TRUE)
    cbind(rep(totdist,101), rep(meanvelo,101))
    })
    data$totdist <- sumdist[,2]
    data$meanvelo <- sumdist[,3]

  
  
  #### step 10 - M: bigflips ####

  f_bigflip <- function(x) {
    if(max(x[,7]) > abs(box.cor$start[1]-box.cor$left[1]))
    { bflip <- 1 } else { bflip <- 0}
    cbind(rep(bflip,101))
  }


  v_bflip <- ddply(data, c("id"), f_bigflip)[,2]
  data$bflip <- v_bflip

  return(data)
}




# call function

dat_processed <- mt.preprocess(d_test, box.cor)





