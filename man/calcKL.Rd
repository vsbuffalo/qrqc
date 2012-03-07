\name{calcKL-methods}
\alias{calcKL}
\alias{calcKL,SequenceSummary-method}
\title{Calculate the Kullback-Leibler Divergence Between the k-mer
  Distribution by Position and the k-mer Distribution Across All
  Positions.
}
\description{
  
  \code{calcKL} takes in an object that inherits from
  \code{SequenceSummary} that has a kmers slot, and returns the terms of
  the K-L divergence sum (which correspond to items in the sample space,
  in this case, k-mers). 

}

\usage{
  calcKL(x)
}

\arguments{
  \item{x}{an S4 object a class that inherits from \code{SequenceSummary}.}
}

\value{
  \code{calcKL} returns a \code{data.frame} with columns:

  \item{kmer}{the k-mer sequence.}
  \item{position}{the position in the read.}
  \item{kl}{the K-L term for this k-mer in the K-L sum, calculated as
    p(i)*log2(p(i)/q(i)).}
  \item{p}{the probability for this k-mer, at this position.}
  \item{q}{the probability for this k-mer across all positions.}
}

\note{
  The K-L divergence calculation in \code{calcKL} uses base 2 in the
  log; the units are in bits.
}

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
  ## Load a somewhat contaminated FASTQ file
  s.fastq <- readSeqFile(system.file('extdata', 'test.fastq',
    package='qrqc'), hash.prop=1)

  ## As with getQual, this function is provided so custom graphics can
  ## be made easily. For example K-L divergence by position:
  kld <- with(calcKL(s.fastq), aggregate(kl, list(position),
    sum))
  colnames(kld) <- c("position", "KL")
  p <- ggplot(kld) + geom_line(aes(x=position, y=KL), color="blue")
  p + scale_y_continuous("K-L divergence")
}

\keyword{methods}
\keyword{graphics}
