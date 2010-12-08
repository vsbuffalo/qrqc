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
  \item{obj}{an S4 object of class \code{SequenceSummary} from
    \code{readSeqFile}.}
}
\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
  ## Load a FASTQ file, with sequence hashing.
  s.fastq <- readSeqFile(system.file('data/test.fastq', package='qrqc'))

  ## Plot Qualities
  plotSeqLengths(s.fastq)
}
\keyword{graphics}