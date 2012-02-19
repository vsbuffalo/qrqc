\name{qualPlot-methods}
\docType{methods}
\alias{qualPlot-methods}
\alias{qualPlot,FASTQSummary-method}
\alias{qualPlot,list-method}
\title{Plot a Base Quality Boxplot by Position}
\description{
  
  \code{basePlot} plots quality statistics by position. Optionally, it
  adds a smoothed curve through the qualities, which is fit with data
  randomly drawn from the distribution of qualities at each position.
  
  The grey lines indicate the range of the 10\% and 90\% quantiles, the
  orange lines indicate the range of the 25\% and 75\% quartiles, the
  blue point is the mean. Optionally, one can plot the median as well.
  
}

\usage{
  qualPlot(x, smooth=TRUE, extreme.color="grey", quartile.color="orange", mean.color="blue", median.color=NULL)
}

\arguments{
  \item{x}{an S4 object of class \code{FASTQSummary} from
    \code{readSeqFile} or a list of \code{FASTQSummary} objects with
    names.}
  \item{smooth}{a logical value indicating whether to add a smooth
    curve.}
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


\section{Methods}{
\describe{

\item{\code{signature(x = "FASTQSummary")}}{
  \code{qualPlot} will plot the qualities for a single object of class \code{FASTQSummary}.
}

\item{\code{signature(x = "list")}}{
  \code{qualPlot} will plot the qualities for each of the
  \code{FASTQSummary} items in the list and display them in a series of
  panels.}
}}

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
  ## Load a FASTQ file, with sequence hashing.
  s.fastq <- readSeqFile(system.file('extdata', 'test.fastq', package='qrqc'))

  ## Plot qualities
  qualPlot(s.fastq)

  ## Combine with ggplot2 to produce custom graphics
  p <- qualPlot(s.fastq, smooth=FALSE)
  p <- p + geom_smooth(aes(x=position, y=quality),
    data=getMCQual(s.fastq), method="lm", color="green", se=FALSE)
  p <- p + theme_bw()
  p

  ## Plot multiple quality plots
  s.trimmed.fastq <- readSeqFile(system.file('extdata',
    'test-trimmed.fastq', package='qrqc'))
  qualPlot(list("not trimmed"=s.fastq, "trimmed"=s.trimmed.fastq))
}

\seealso{getQual}
\keyword{methods}
\keyword{graphics}
