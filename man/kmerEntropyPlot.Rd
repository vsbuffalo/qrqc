\name{kmerEntropyPlot-methods}
\docType{methods}
\alias{kmerEntropyPlot}
\alias{kmerEntropyPlot-methods}
\alias{kmerEntropyPlot,SequenceSummary-method}
\alias{kmerEntropyPlot,list-method}
\title{Plot Entropy of k-mers by Position}
\description{
  
  \code{kmerEntropyPlot} plots the Shannon entropy per position of
  k-mers. Lower Shannon entropy implies that the distribution of k-mers
  is non-random and could indicate bias.
}

\usage{
  kmerEntropyPlot(x)
}

\arguments{
  \item{x}{an S4 object a class that inherits from \code{SequenceSummary} from
    \code{readSeqFile} or a list of objects that inherit from \code{SequenceSummary} with
    names.}
}


\section{Methods}{
\describe{
  \item{\code{signature(x = "SequenceSummary")}}{
    \code{kmerEntropyPlot} will plot Shannon entropy per position for an
    object that inherits from \code{SequenceSummary}.
  }
  
  \item{\code{signature(x = "list")}}{
    \code{kmerEntropyPlot} will plot the Shannon entropy per position
    for each of the objects that inherit from \code{SequenceSummary} in
    the list and display them in a series of panels.
  }
}}

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
  ## Load a somewhat contaminated FASTQ file
  s.fastq <- readSeqFile(system.file('extdata', 'test.fastq',
    package='qrqc'), hash.prop=1)

  ## Load a really contaminated FASTQ file
  s.contam.fastq <- readSeqFile(system.file('extdata',
    'test-contam.fastq', package='qrqc'), hash.prop=1)

  ## Load a random (equal base frequency) FASTA file
  s.random.fasta <- readSeqFile(system.file('extdata',
    'random.fasta', package='qrqc'), type="fasta", hash.prop=1)

  ## Plot the Shannon entropy for a single file
  kmerEntropyPlot(s.fastq)

  ## Plot the Shannon entropy for many files
  kmerEntropyPlot(list("highly contaminated"=s.contam.fastq, "less
    contaminated"=s.fastq, "random"=s.random.fasta))
}

\seealso{\code{\link{getKmer}}, \code{\link{calcKL}}, \code{\link{kmerKLPlot}}}
\keyword{methods}
\keyword{graphics}
