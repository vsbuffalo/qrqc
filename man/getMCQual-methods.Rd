\name{getMCQual-methods}
\docType{methods}
\alias{getMCQual}
\alias{getMCQual-methods}
\alias{getMCQual,FASTQSummary-method}
\title{Get a Data Frame of Simulated Qualitied from a \code{FASTQSummary} object}
\description{
  An object that inherits from class \code{\link[=FASTQSummary-class]{FASTQSummary}} contains
  base quality data by position gathered by
  \code{\link{readSeqFile}}. \code{\link{getMCQual}} generates simulated
  quality data for each base from this binned  quality data that can be
  used for adding smoothed lines via lowess.  

  This accessor function is useful if you want to map variables to
  custom \code{ggplot2} aesthetics.
}


\usage{
  getMCQual(x, n=100)
}

\arguments{
  \item{x}{an S4 object that inherits from \code{FASTQSummary} from
    \code{readSeqFile}.}
  \item{n}{a numeric value indicating the number of quality values to
    draw per base.}
}


\value{
  \code{getMCQual} returns a \code{data.frame} with columns:

  \item{position}{the position in the read.}
  \item{quality}{simulated quality scores.}
}

\section{Methods}{
  \describe{
    \item{\code{signature(x = "FASTQSummary")}}{
      \code{getMCQual} is a function that works on any object with class
      \code{FASTQSummary} read in with \code{readSeqFile}.
    }
}}

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
  ## Load a FASTQ file, with sequence hashing.
  s.fastq <- readSeqFile(system.file('extdata', 'test.fastq',
    package='qrqc'))

  # A custom quality plot
  ggplot(getQual(s.fastq)) + geom_linerange(aes(x=position, ymin=lower,
    ymax=upper), color="grey") + geom_smooth(aes(x=position, y=quality),
    data=getMCQual(s.fastq), color="blue", se=FALSE)
}
\seealso{\code{\link{getGC}}, \code{\link{getSeqlen}}, \code{\link{getBase}}, \code{\link{getBaseProp}}, \code{\link{getQual}}, \code{\link{qualPlot}}}

\keyword{methods}
\keyword{accessor}
