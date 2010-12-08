\name{plotBaseProps}
\alias{plotBaseProps}
\title{Plot Proportions of Bases by Position}
\description{
  \code{plotBaseProps} plots the proporition of bases by position
  (cycle) in the read.
}
\usage{
  plotBaseProps(obj)
}
\arguments{
  \item{obj}{an S4 object of class \code{SequenceSummary} from \code{readSeqFile}.}
}
\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
  ## Load a FASTQ file, with sequence hashing.
  s.fastq <- readSeqFile(system.file('data/test.fastq', package='qrqc'))

  ## Plot Qualities
  plotBaseProps(s.fastq)
}
\keyword{graphics}