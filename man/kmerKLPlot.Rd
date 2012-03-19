\name{kmerKLPlot-methods}
\docType{methods}
\alias{kmerKLPlot}
\alias{kmerKLPlot-methods}
\alias{kmerKLPlot,SequenceSummary-method}
\alias{kmerKLPlot,list-method}
\title{Plot K-L Divergence Components for a Subset of k-mers to Inspect for Contamination}
\description{
  
  \code{kmerKLPlot} calls \code{\link{calcKL}}, which calculates the
  Kullback-Leibler divergence between the k-mer distribution at each
  position compared to the k-mer distribution across all
  positions. \code{kmerKLPlot} then plots each k-mer's contribution to
  the total K-L divergence by stack bars, for a \emph{subset} of the
  k-mers. Since there are 4^k possible k-mers for some value k-mers,
  plotting each often dilutes the interpretation; however one can
  increase \code{n.kmers} to a number greater than the possible number
  of k-mers to force \code{kmerKLPlot} to plot the entire K-L divergence
  and all terms (which are k-mers) in the sum.

  If a \code{x} is a \code{list}, the K-L k-mer plots are faceted by
  sample; this allows comparison to a FASTA file of random reads.

  Again, please note that this is \emph{not} the total K-L divergence,
  but rather the K-L divergence calculated on a subset of the sample
  space (those of the top \code{n.kmers} k-mers selected).
}

\usage{
  kmerKLPlot(x, n.kmers=20)
}

\arguments{
  \item{x}{an S4 object a class that inherits from \code{SequenceSummary} from
    \code{readSeqFile} or a list of objects that inherit from \code{SequenceSummary} with
    names.}
  \item{n.kmers}{a integer value indicating the size of top k-mers to include.}
}


\section{Methods}{
\describe{
  \item{\code{signature(x = "SequenceSummary")}}{
    \code{kmerKLPlot} will plot the K-L divergence for a subset of k-mers for a single object that
    inherits from \code{SequenceSummary}.
  }
  
  \item{\code{signature(x = "list")}}{
    \code{kmerKLPlot} will plot the K-L divergence for a susbet of
    k-mers  for each of the objects that inherit from 
    \code{SequenceSummary} in the list and display them in a series of
    panels.}
}}

\note{

  The K-L divergence calculation in \code{calcKL} uses base 2 in the
  log; the units are in bits.

  Also, note that \code{ggplot2} warns that "Stacking is not well defined when ymin
  != 0". This occurs when some k-mers are less frequent in the positional
  distribution than the distribution across all positions, and the term of
  the K-L sum is negative (producing a bar below zero). This does not
  appear to affect the plot much. In examples below, warnings are
  suppressed, but the given this is a valid concern from \code{ggplot2},
  warnings are not suppressed in the function itself.
  
}

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

  ## Make K-L divergence plot - shows slight 5'-end bias. Note units
  ## (bits)
  suppressWarnings(kmerKLPlot(s.fastq))

  ## Plot multiple K-L divergence plots
  suppressWarnings(kmerKLPlot(list("highly contaminated"=s.contam.fastq, "less
    contaminated"=s.fastq, "random"=s.random.fasta)))
}
\seealso{\code{\link{getKmer}}, \code{\link{calcKL}},
  \code{\link{kmerEntropyPlot}}}
\keyword{methods}
\keyword{graphics}
