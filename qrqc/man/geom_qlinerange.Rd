\name{geom_qlinerange}
\alias{geom_qlinerange}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
geom_qlinerange(extreme.color = "grey", quantile.color = "orange", mean.color = "blue", median.color = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{extreme.color}{
%%     ~~Describe \code{extreme.color} here~~
}
  \item{quantile.color}{
%%     ~~Describe \code{quantile.color} here~~
}
  \item{mean.color}{
%%     ~~Describe \code{mean.color} here~~
}
  \item{median.color}{
%%     ~~Describe \code{median.color} here~~
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
function (extreme.color = "grey", quantile.color = "orange", 
    mean.color = "blue", median.color = NULL) 
{
    args <- as.list(match.call(call = sys.call(sys.parent()))[-1])
    l <- list(extreme.color = geom_linerange(aes(x = position, 
        ymin = alt.lower, ymax = alt.upper), color = extreme.color), 
        quantile.color = geom_linerange(aes(x = position, ymin = lower, 
            ymax = upper, y = mean), color = quantile.color, 
            size = 1.2), mean.color = geom_point(aes(x = position, 
            y = mean), color = mean.color), median.color = geom_point(aes(x = position, 
            y = middle), color = median.color))
    keep <- sapply(names(l), function(n) !is.null(get(n)))
    l[keep]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
