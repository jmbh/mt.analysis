####################################################
### Preprocess mouse-tracking data for analysis
### jonashaslbeck@gmail.com
### January 2015
####################################################

#depedencies
#library(plyr)

## input:
# 1) data matrix containing the columns: trialid defining each trajectory uniquely,x,y,time,chosen box (left/right)
# 2) coordinates of the starting position and the two boxes of the mouse


## output:
# 1) time-normalized x and y values lined up to (0,0); x also fliped on one side
# 2) on time step level: velocity(difference), distance to direct line
# 3) on trial level: RT, mean velocity, MAD, AAD, directionalchange(DC), total distance

##TESTINPUT
#box.cor <- list("start"=c(960,230), "left"=c(130,905), "right"=c(1830,905)) #coordinates
#i.id <- c("part_id", "trial_count", "sample_count")
#i.measure <- c("x", "y","time", "box_chosen")


mt.preprocess <- function(
  data, # n x p data matrix
  box.cor, # coordinates of start position and boxes, see details
  i.id, # names of columns whose combination indicate unique trajectories
  i.measure, #names of columns indicating x,y,time,chosen box in this order
  tsteps = 101 #numer of timesteps, default = 101
  ) {
  
  #### step 1 - sanity checks 

  if(length(box.cor) != 3) {
    stop("Please specify the positions of the boxes.")
    }
  if(4!=(length(i.measure))) {
    stop("Pleas specify i.measure correctly.")
    }
  if(sum(is.na(data)>0)) {
    stop("No missing values allowed.")
    }
  
  num_check <- apply(data[,i.measure], 2, is.numeric) == FALSE
  
  if(sum(num_check)>0) {
    stop("Only numerical values allowed.")
  }
  

data <- as.data.frame(data)
ids <- length(i.id)
cn <- colnames(data)

#reorder data

data <- cbind(data[,i.measure[1]], data[,i.measure[2]], data[,i.measure[3]], data[,i.measure[4]], 
              data[,(cn %in% i.measure)==FALSE])
colnames(data)[1:4] <- c("x", "y", "t", "b")
data$x <- as.numeric(data$x)
data$y <- as.numeric(data$y)
data$t <- as.numeric(data$t)

#### step 2 -  set starting point to zero  

  f_settozero <- function(z) 
  {
    x <- z[1] - box.cor$start[1]
    y <- z[2] - box.cor$start[2]
    return(rbind(x,y))
  }
  
  new.cord <- t(apply(data[,c("x", "y")], 1, f_settozero))
  data$x <- new.cord[,1]
  data$y <- new.cord[,2]
  


  #### step 3 - calc reaction time

  f_calc_rts <- function(x) 
  {
    le <- length(x$t)
    rt <- x$t[le] - x$t[1]
    rt.var <- cbind(rep(rt, le))
    return(rt.var)
  }

  data$rt <- ddply(data, c(i.id), f_calc_rts)[,(ids+1)]


#### step 4 - time normalise ####

  f_timenorm <- function(x) 
  {
    v_x <- x[,1]
    v_y <- x[,2]
    v_time <- x[,3]
    
    #normalise time vector [0,1]
    v_time_norm <- (v_time - v_time[1]) / (v_time[length(v_time)]-v_time[1])
    v_time_norm <- v_time_norm*(tsteps-1)
    
    # interpolation of x&y to tsteps equally spaced time slices
    lin.x <- approx(v_time_norm, v_x, xout = 0:(tsteps-1), method = "linear") # interpolate x coordinates
    lin.y <- approx(v_time_norm, v_y, xout = 0:(tsteps-1), method = "linear") # interpolate y coordinates
    
    rawdata_restoftable <- x[rep(1,tsteps),4:ncol(data)]

    
    data_export <- cbind(lin.x$y, lin.y$y, lin.x$x, rawdata_restoftable) 
    return(data_export)
  }

  data_norm <- ddply(data, c(i.id), f_timenorm)
  colnames(data_norm)[1:3] <- c("x", "y", "t")
  data <- data_norm

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
  
  xflip <- apply(cbind(data[,c("b")], data[,c("x")]),1,f_flip)
  data$xflip <- xflip


  #### step 6 - M: distance to direct line ####

  f_calculatedistance<- function(z)
  {
          
    slopedl <- (box.cor$left[2]-box.cor$start[2])/(box.cor$left[1]-box.cor$start[1]) #slope of direct line
    
    #find line with slope = -1/slopedl, which passes through the empirical point (interesting thing is only the intercept with y-axis, because slopes are known)
    c_x<-as.numeric(z["xflip"])
    c_y<-as.numeric(z["y"])
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
    AAD <- rep(mean(x$dist),tsteps) # aad
    MAD <- rep(max(x$dist),tsteps) # mad
    out <- cbind(AAD, MAD)
    return(out)
  }
  
  v_aadmad <- ddply(data, c(i.id), f_aadmad)
  data$AAD <- v_aadmad$AAD
  data$MAD <- v_aadmad$MAD
  

  #### step 8 - M: Velocity ####
  
  dx <- c(0,data$xflip[-length(data$xflip)]) #dummy variables, vector shifted by one
  dy <- c(0,data$y[-length(data$y)])
  velo <- sqrt((data$xflip-dx)^2+(data$y-dy)^2) #euclidean distance between (x,y) at t-1 and t
  
  
  data$velo <- velo
  data$velo[data$t==0] <- NA #there cant be velocity in the first timestamp (no reference before)
  

  #### step 9 - M: Total distance / mean Velocity ####
    sumdist <- ddply(data, c(i.id), function(x) {
    totdist <- sum(x$velo, na.rm=TRUE)
    meanvelo <- mean(x$velo, na.rm=TRUE)
    cbind(rep(totdist,tsteps), rep(meanvelo,tsteps))
    })
    data$totdist <- sumdist[,(ncol(sumdist)-1)]
    data$meanvelo <- sumdist[,ncol(sumdist)]

  

  #re-order dataframe
  data_out <- cbind(data[,c(i.id)], data[,colnames(data)[(colnames(data) %in% i.id)!=TRUE]])

  return(data_out)
} #end of function


## testing function
#t1 <- proc.time()[1]
#dat_processed <- mt.preprocess(dat1, box.cor, i.id, i.measure)
#proc.time()[1] - t1

getwd()
load("data/dataraw.RData")

head(dataraw) #example dataset
box.cor <- list("start"=c(960,230), "left"=c(130,905), "right"=c(1830,905)) #coordinates from the example dataset
i.id <- c("trial")
i.measure <- c("x", "y","t", "b")
data.norm <- mt.preprocess(data=dataraw,
                      box.cor = box.cor, 
                      i.id, 
                      i.measure, 
                      tsteps=101)
head(data.norm[,1:9])


