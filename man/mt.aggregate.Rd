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


head(dataprocessed) #example dataset

data.aggregated <- mt.aggregate(data = dataprocessed, 
                                i.aggr=c("group"), 
                                i.xyt=c("x", "y", "t"))
                                
head(data.aggregated)                    

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
