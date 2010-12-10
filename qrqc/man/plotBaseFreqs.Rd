\name{plotBaseFreqs}
\alias{plotBaseFreqs}
\alias{plotBaseFreqs,SequenceSummary-method}
\alias{plotBaseFreqs,SequenceSummary,character-method}
\alias{plotBaseFreqs,SequenceSummary,missing-method}
\title{Plot Frequencies of Bases by Position}
\description{
  \code{plotBaseFreqs} plots the proportion of bases by position in the
  read.
}
\usage{
  plotBaseFreqs(obj, bases=NULL)
}
\arguments{
  \item{obj}{an S4 object of class that inherits from \code{\linkS4class{SequenceSummary}} (either 
    \code{\linkS4class{FASTASummary}} or \code{\linkS4class{FASTQSummary}}) from
    \code{readSeqFile}.}
  \item{bases}{a vector of characters indicating which bases to
    include. The default value \code{NULL} indicates to plot _all_ bases.}
}
\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
  ## Load a FASTQ file, with sequence hashing.
  s.fastq <- readSeqFile(system.file('extdata/test.fastq', package='qrqc'))

  ## Plot Qualities
  plotBaseFreqs(s.fastq)
}
\keyword{graphics}