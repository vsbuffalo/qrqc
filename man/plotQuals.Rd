\name{plotQuals-methods}
\docType{methods}
\alias{plotQuals}
\alias{plotQuals-methods}

\title{Plot a Base Quality Boxplot by Position}
\description{
  \code{plotQuals} plots quality statistics by position.
  Optionally, it adds a lowess curve through the qualities,
  which is fit with data randomly drawn from the distribution of
  qualities at each position. A histogram of the sequence length
  distribution is plotted above the quality plot when \code{histogram}
  is \code{TRUE}.

  A legend is plotted on the bottom left if \code{legend} is
  \code{TRUE} (this location is used because this where the bases are
  likely to be of highest quality, and thus not overlap the legend). The
  grey lines indicate the range of the 10\% and 90\% quantiles, the
  orange lines indicate the range of the 25\% and 75\% quartiles, the
  blue point is the median, the green dash is the mean, and the purple
  line is the lowess curve if \code{lowess} is \code{TRUE}.
  
}
\usage{
  plotQuals(obj, ylim='relative', lowess=TRUE, histogram=TRUE, legend=TRUE)
}
\arguments{
  \item{obj}{an S4 object of class \code{FASTQSummary} from \code{readSeqFile}.}
  \item{ylim}{either `relative' or `fixed', which will scale the y axis
    to either the relative range (from the data) or absolute range of qualities.}
  \item{lowess}{a logical value indicating whether to fit a lowess curve
    through the quality plot.}
  \item{histogram}{a logical value indicating whether to add a histogram
    of the sequence length distribution above the quality plot.}
  \item{legend}{a logical value indicating whether a legend is to be included.}
}
\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
\dontrun{
  ## Load a FASTQ file, with sequence hashing.
  s.fastq <- readSeqFile(system.file('extdata', 'test.fastq', package='qrqc'))

  ## Plot Qualities
  plotQuals(s.fastq)
}
}
\seealso{qualPlot}
\keyword{graphics}
\keyword{methods}
