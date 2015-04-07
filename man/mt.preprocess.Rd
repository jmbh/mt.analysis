\name{mt.preprocess}
\alias{mt.preprocess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Peprocess motor-tracing data
}
\description{
Time-normalizes the trajectories and computes a set of derivatives of the trajectories}
}
\usage{
mt.preprocess(data, box.cor, i.id, i.measure, ts = 101)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
n x m - matrix with n measurements and m variables
}
  \item{box.cor}{
%%     ~~Describe \code{box.cor} here~~
}
  \item{i.id}{
%%     ~~Describe \code{i.id} here~~
}
  \item{i.measure}{
%%     ~~Describe \code{i.measure} here~~
}
  \item{ts}{
%%     ~~Describe \code{ts} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, box.cor, i.id, i.measure, ts = 101) 
{
    if (length(box.cor) != 3) {
        stop("Please specify the positions of the boxes.")
    }
    if (4 != (length(i.measure))) {
        stop("Pleas specify i.measure correctly.")
    }
    if (sum(is.na(data) > 0)) {
        stop("No missing values allowed.")
    }
    data <- as.data.frame(data)
    ids <- length(i.id)
    cn <- colnames(data)
    data <- cbind(data[, i.measure[1]], data[, i.measure[2]], 
        data[, i.measure[3]], data[, i.measure[4]], data[, (cn \%in\% 
            i.measure) == FALSE])
    colnames(data)[1:4] <- c("x", "y", "t", "b")
    data$x <- as.numeric(data$x)
    data$y <- as.numeric(data$y)
    data$t <- as.numeric(data$t)
    f_settozero <- function(z) {
        x <- z[1] - box.cor$start[1]
        y <- z[2] - box.cor$start[2]
        return(rbind(x, y))
    }
    new.cord <- t(apply(data[, c("x", "y")], 1, f_settozero))
    data$x <- new.cord[, 1]
    data$y <- new.cord[, 2]
    f_calc_rts <- function(x) {
        le <- length(x$t)
        rt <- x$t[le] - x$t[1]
        rt.var <- cbind(rep(rt, le))
        return(rt.var)
    }
    data$rt <- ddply(data, c(i.id), f_calc_rts)[, (ids + 1)]
    hist(data$rt)
    f_timenorm <- function(x) {
        v_x <- x[, 1]
        v_y <- x[, 2]
        v_time <- x[, 3]
        v_time_norm <- (v_time - v_time[1])/(v_time[length(v_time)] - 
            v_time[1])
        v_time_norm <- v_time_norm * (ts - 1)
        lin.x <- approx(v_time_norm, v_x, xout = 0:(ts - 1), 
            method = "linear")
        lin.y <- approx(v_time_norm, v_y, xout = 0:(ts - 1), 
            method = "linear")
        rawdata_restoftable <- x[rep(1, 101), 4:ncol(data)]
        data_export <- cbind(lin.x$y, lin.y$y, lin.x$x, rawdata_restoftable)
        return(data_export)
    }
    data_norm <- ddply(data, c(i.id), f_timenorm)
    colnames(data_norm)[1:3] <- c("x", "y", "t")
    data <- data_norm
    f_flip <- function(z) {
        x <- z[2]
        if (z[1] == 1) {
            x <- x * (-1)
        }
        return(x)
    }
    xflip <- apply(cbind(data[, c("b")], data[, c("x")]), 1, 
        f_flip)
    data$xflip <- xflip
    f_calculatedistance <- function(z) {
        slopedl <- (box.cor$left[2] - box.cor$start[2])/(box.cor$left[1] - 
            box.cor$start[1])
        c_x <- as.numeric(z["xflip"])
        c_y <- as.numeric(z["y"])
        slopedl_perp <- -1/slopedl
        c <- -slopedl_perp * c_x + c_y
        x_intersect <- (-0 + c)/(slopedl - slopedl_perp)
        y_intersect <- x_intersect * slopedl
        c_distance <- sqrt((x_intersect - c_x)^2 + (y_intersect - 
            c_y)^2)
        return(c_distance)
    }
    data$dist <- apply(data, 1, f_calculatedistance)
    f_aadmad <- function(x) {
        AAD <- rep(mean(x$dist), 101)
        MAD <- rep(max(x$dist), 101)
        out <- cbind(AAD, MAD)
        return(out)
    }
    v_aadmad <- ddply(data, c(i.id), f_aadmad)
    data$AAD <- v_aadmad$AAD
    data$MAD <- v_aadmad$MAD
    dx <- c(0, data$xflip[-length(data$xflip)])
    dy <- c(0, data$y[-length(data$y)])
    velo <- sqrt((data$xflip - dx)^2 + (data$y - dy)^2)
    data$velo <- velo
    data$velo[data$t == 0] <- NA
    sumdist <- ddply(data, c(i.id), function(x) {
        totdist <- sum(x$velo, na.rm = TRUE)
        meanvelo <- mean(x$velo, na.rm = TRUE)
        cbind(rep(totdist, 101), rep(meanvelo, 101))
    })
    data$totdist <- sumdist[, (ncol(sumdist) - 1)]
    data$meanvelo <- sumdist[, ncol(sumdist)]
    data_out <- cbind(data[, c(i.id)], data[, colnames(data)[(colnames(data) \%in\% 
        i.id) != TRUE]])
    return(data_out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
