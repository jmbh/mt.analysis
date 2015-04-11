\name{mt.aggregate}
\alias{mt.aggregate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Aggregate trajectories
}
\description{
Aggregates trajectories according to a factorial design based on specified variables.}
\usage{
mt.aggregate(data, i.aggr, i.xyt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
n x m - matrix with n measurements and m variables
}
  \item{i.aggr}{
Vector of column names of the variables used for aggregation (if there is more than one variable the data is aggregated for each cell of the factorial design)
}
  \item{i.xyt}{
Column names of columns containing x,y and time (in that order)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a matrix with ....
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


> head(dataprocessed) #example dataset
  trial       x       y t group  rt   xflip     dist      AAD      MAD      velo  totdist meanvelo
1     1 127.000 -44.000 0     1 466 127.000 45.99351 73.89085 109.7982        NA 1191.295 11.91295
2     1 125.536 -41.804 1     1 466 125.536 46.77353 73.89085 109.7982 2.6392635 1191.295 11.91295
3     1 125.000 -41.000 2     1 466 125.000 47.05911 73.89085 109.7982 0.9662877 1191.295 11.91295
4     1 123.408 -40.204 3     1 466 123.408 46.67220 73.89085 109.7982 1.7799101 1191.295 11.91295
5     1 123.000 -40.000 4     1 466 123.000 46.57305 73.89085 109.7982 0.4561579 1191.295 11.91295
6     1 123.000 -40.000 5     1 466 123.000 46.57305 73.89085 109.7982 0.0000000 1191.295 11.91295


data.aggregated <- mt.aggregate(data = dataprocessed, 
                                i.aggr=c("group"), 
                                i.xyt=c("x", "y", "t"))
                                
> head(data.aggregated)                    
  group t        mx        my
1     0 0 110.35000 -53.75000
2     0 1 108.69845 -51.59923
3     0 2 106.04540 -48.80920
4     0 3 103.74498 -46.29978
5     0 4  99.03794 -41.29537
6     0 5  93.54546 -35.21787

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
