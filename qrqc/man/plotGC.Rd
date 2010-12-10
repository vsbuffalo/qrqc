\name{plotGC}
\alias{plotGC}
\title{Plot per Base GC Content by Position}
\description{
  \code{plotGC} plots the GC proportion by position.
}
\usage{
  plotGC(obj)
}
\arguments{
  \item{obj}{an S4 object of class that inherits from \code{\linkS4class{SequenceSummary}} (either 
    \code{\linkS4class{FASTASummary}} or \code{\linkS4class{FASTQSummary}}) from
    \code{readSeqFile}.}
}
\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
  ## Load a FASTQ file, with sequence hashing.
  s.fastq <- readSeqFile(system.file('extdata/test.fastq', package='qrqc'))

  ## Plot Qualities
  plotGC(s.fastq)
}
\keyword{graphics}