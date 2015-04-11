\name{mt.ttdetect}
\alias{mt.ttdetect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Detect different trajectory types
}
\description{
Uses k-means clustering on a specified vector of variables to identify different types of trajectories
}
\usage{
mt.ttdetect(data, kclust, varclust)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The dataset containing the mouse-tracing data and the variables used for clustering
}
  \item{kclust}{
Number of clusters
}
  \item{varclust}{
Vector of variable names used for clustering
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

head(dataprocessed) #example dataset


data.g <- at.ttdetect(dataprocessed, 
                      kclust = 2, 
                      varclust = "MAD")
                      
head(data.g)
table(data.g$clust)





}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
