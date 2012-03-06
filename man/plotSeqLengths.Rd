\name{plotSeqLengths-methods}
\docType{methods}
\alias{plotSeqLengths}
\alias{plotSeqLengths-methods}
\title{Plot Histogram of Sequence Lengths}
\description{
  \code{plotSeqLengths} plots histogram of sequence lengths.
}
\usage{
  plotSeqLengths(obj)
}
\arguments{
  \item{obj}{an S4 object of class that inherits from \code{\linkS4class{SequenceSummary}} (either 
    \code{\linkS4class{FASTASummary}} or \linkS4class{FASTQSummary}) from
    \code{readSeqFile}.}
}
\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
\dontrun{
  ## Load a FASTQ file, with sequence hashing.
  s.fastq <- readSeqFile(system.file('extdata', 'test.fastq', package='qrqc'))

  ## Plot Qualities
  plotSeqLengths(s.fastq)
}
}
\seealso{seqlenPlot}
\keyword{graphics}