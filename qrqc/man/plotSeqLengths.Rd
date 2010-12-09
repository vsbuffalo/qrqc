\name{plotSeqLengths}
\alias{plotSeqLengths}
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
  ## Load a FASTQ file, with sequence hashing.
  s.fastq <- readSeqFile(system.file('extdata/test.fastq', package='qrqc'))

  ## Plot Qualities
  plotSeqLengths(s.fastq)
}
\keyword{graphics}