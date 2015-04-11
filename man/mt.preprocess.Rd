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
mt.preprocess(data, box.cor, i.id, i.measure, tsteps = 101)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
n x m - matrix with n measurements and m variables
}
  \item{box.cor}{
List containting the position of the center of the two response boxes and the starting position of the cursor
}
  \item{i.id}{
Vector of column names indicating the id-variables identifying unique trajectries
}
  \item{i.measure}{
Vector of column names indicating the columns for the measures of x, y, time and the selected box (left=0 or right=1)
}
  \item{tsteps}{
The number of time-steps the data is normalized to. Default = 101.
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
Spivey, M. J., Grosjean, M., & Knoblich, G. (2005). Continuous attraction toward phonological competitors. Proceedings of the National Academy of Sciences of the United States of America, 102(29), 10393-10398.
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
