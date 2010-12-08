\name{plotQuals}
\alias{plotQuals}
\title{Plot a Base Quality Boxplot by Position}
\description{
  \code{plotQuals} plots quality statistics by position
  (cycle). Optionally, it adds a lowess curve through the qualities,
  which is fit with data randomly drawn from the distribution of
  qualities at each position. A histogram of the sequence length
  distribution is plotted above the quality plot when \code{histogram}
  is \code{TRUE}.
}
\usage{
  plotQuals(obj, ylim='relative', lowess=TRUE, histogram=TRUE)
}
\arguments{
  \item{obj}{an S4 object of class \code{SequenceSummary} from \code{readSeqFile}.}
  \item{ylim}{either `relative' or `fixed', which will scale the y axis
    to either the relative range (from the data) or absolute range of qualities.}
  \item{lowess}{a logical value indicating whether to fit a lowess curve
    through the quality plot.}
  \item{histogram}{a logical value indicating whether to add a histogram
  of the sequence length distribution above the quality plot.}
}
\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
  ## Load a FASTQ file, with sequence hashing.
  s.fastq <- readSeqFile(system.file('data/test.fastq', package='qrqc'))

  ## Plot Qualities
  plotQuals(s.fastq)
}
\keyword{graphics}