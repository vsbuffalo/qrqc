\name{plotBases}
\alias{plotBases}
\title{Plot Bases by Position}
\description{
  
  \code{plotBases} plots the frequency or proportion of bases by
  position in the read.

  \code{plotBases} uses the Sanger base color scheme: blue is Cytosine,
  green is Adenine, black is Guanine, red is Thymine, and purple in N
  (any base). When further IUPAC nucleotides are found in the sequences,
  a legend is added.
  
}
\usage{
  plotBases(obj, bases=NULL, type="freq")
}
\arguments{
  
  \item{obj}{an S4 object of class that inherits from
    \code{\linkS4class{SequenceSummary}} (either
    \code{\linkS4class{FASTASummary}} or
    \code{\linkS4class{FASTQSummary}}) from \code{readSeqFile}.}
  
  \item{bases}{a vector of characters indicating which bases to
    include. The default value \code{NULL} indicates to plot _all_
    bases.}
  
  \item{type}{a character string that is either "freq" or "prop"
    indicating whether to plot frequencies or proportions on the
    y-axis.}

}
\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
  ## Load a FASTQ file, with sequence hashing.
  s.fastq <- readSeqFile(system.file('extdata', 'test.fastq', package='qrqc'))

  ## Plot base frequencies
  plotBases(s.fastq, type="freq")

  ## Plot base proportions
  plotBases(s.fastq, type="prop")
}
\keyword{graphics}