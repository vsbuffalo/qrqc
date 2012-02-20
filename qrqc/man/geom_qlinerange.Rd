\name{geom_qlinerange}
\alias{geom_qlinerange}
\title{Use Line Segments and Points to Plot Quality Statistics by Position in the Read}
\description{
  \code{geom_qlinerange} uses multiple line segments and points to plot
  quality ranges. By default the 10\% and 90\% range in plotted in grey, the quartile
  range in orange, and the mean as a point in blue.
}

\usage{
geom_qlinerange(extreme.color="grey", quartile.color="orange", mean.color="blue", median.color=NULL)
}

\arguments{
  \item{extreme.color}{a character value indicating the color to use for
    the extreme values (the 10\% and 90\% quantiles). If \code{NULL},
    these line seqments will not be added.}
  \item{quartile.color}{a character value indicating the color to use
    for the quartiles.  If \code{NULL}, these line seqments will not be
    added.}
  \item{mean.color}{a character value indicating the color to use for
    the mean. If \code{NULL}, these line seqments will not be added.}
  \item{median.color}{a character value indicating the color to use for
    the median. If \code{NULL}, these line seqments will not be added.}
}

\value{
  A list of geoms from \code{ggplot2} that this function put together,
  to be added to a call to \code{ggplot} which contains a data frame of
  quality data, i.e. from \code{getQual}.
}

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\seealso{getQual}
\keyword{methods}
\keyword{graphics}
