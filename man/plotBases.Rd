\name{plotBases-methods}
\docType{methods}
\alias{plotBases}
\alias{plotBases-methods}
\title{Plot Bases by Position}
\description{
  
  \code{plotBases} plots the frequency or proportion of bases by
  position in the read.

  \code{plotBases} uses the Sanger base color scheme: blue is Cytosine,
  green is Adenine, black is Guanine, red is Thymine, and purple in N
  (any base). Other IUPAC nucleotides are colored using \pkg{RColorBrewer}.
  
}
\usage{
  plotBases(obj, type="freq", bases=NULL, legend=TRUE)
}
\arguments{
  
  \item{obj}{an S4 object of class that inherits from
    \code{\linkS4class{SequenceSummary}} (either
    \code{\linkS4class{FASTASummary}} or
    \code{\linkS4class{FASTQSummary}}) from \code{readSeqFile}.}
  
  \item{type}{a character string that is either "freq" or "prop"
    indicating whether to plot frequencies or proportions on the
    y-axis.}
  
  \item{bases}{a vector of characters indicating which bases to
    include. The default value \code{NULL} indicates to plot _all_
    bases.}
  
  \item{legend}{a logical value indicating whether to include a legend on the
  top right.}
  

}
\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
\dontrun{
  ## Load a FASTQ file, with sequence hashing.
  s.fastq <- readSeqFile(system.file('extdata', 'test.fastq', package='qrqc'))

  ## Plot base frequencies
  plotBases(s.fastq, type="freq")

  ## Plot base proportions
  plotBases(s.fastq, type="prop")
}
}
\seealso{basePlot}
\keyword{graphics}