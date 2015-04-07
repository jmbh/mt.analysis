\name{mt.aggregate}
\alias{mt.aggregate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Aggregate trajectories
}
\description{
Aggregates trajectories using a given factor}
\usage{
mt.aggregate(data, i.aggr, i.xyt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
n x m - matrix with n measurements and m variables
}
  \item{i.aggr}{
Column-name of the variable used for aggregation
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, i.aggr, i.xyt) 
{
    out2 <- ddply(data, i.aggr, function(z) {
        out1 <- ddply(z, c("t"), function(x) {
            mx <- mean(as.matrix(x[i.xyt[1]]))
            my <- mean(as.matrix(x[i.xyt[2]]))
            return(cbind(mx, my))
        })
        return(out1)
    })
    return(out2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
