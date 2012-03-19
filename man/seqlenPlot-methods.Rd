\name{seqlenPlot-methods}
\docType{methods}
\alias{seqlenPlot}
\alias{seqlenPlot-methods}
\alias{seqlenPlot,list-method}
\alias{seqlenPlot,SequenceSummary-method}
\title{Plot a Histogram of Sequence Lengths}
\description{
  \code{seqlenPlot} plots a histogram of sequence lengths.
}

\usage{
  seqlenPlot(x)
}

\arguments{
  \item{x}{an S4 object that inherits from \code{SequenceSummary} from
    \code{readSeqFile}.}
}

\section{Methods}{
\describe{
  \item{\code{signature(x = "FASTQSummary")}}{
    \code{seqlenPlot} will plot a histogram of a single object that
    inherits from \code{SequenceSummary}.
  }
  
  \item{\code{signature(x = "list")}}{
    \code{seqlenPlot} will plot a histogram for each of the
    \code{SequenceSummary} items in the list and display them in a
    series of panels.
  }
}}

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
  ## Load a FASTQ file
  s.fastq <- readSeqFile(system.file('extdata', 'test.fastq', package='qrqc'))

  ## Plot sequence lengths
  seqlenPlot(s.fastq)

  ## Plot sequence legnths before and after trimming
  s.trimmed.fastq <- readSeqFile(system.file('extdata',
    'test-trimmed.fastq', package='qrqc'))
  seqlenPlot(list("not trimmed"=s.fastq, "trimmed"=s.trimmed.fastq))
}

\seealso{\code{\link{getSeqlen}}}
\keyword{methods}
\keyword{graphics}

