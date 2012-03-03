\name{getKmer-methods}
\docType{methods}
\alias{getKmer}
\alias{getKmer-methods}
\alias{getKmer,SequenceSummary-method}

\title{Get a Data Frame of k-mer Frequency by Position from a \code{SequenceSummary} Object}

\description{
  An object that inherits from class \code{SequenceSummary} contains
  k-mer frequency data by position gathered by \code{readSeqFile} when
  \code{kmer=TRUE}. \code{getKmer} is an accessor function that is
  useful for custom \code{ggplot2} aesthetics.
}

\usage{
getKmer(x)
}

\arguments{
  \item{x}{an S4 object that inherits from class \code{SequenceSummary}
  from, as returned from \code{readSeqFile}.}
}

\value{
  \code{getKmer} returns a \code{data.frame} with columns:

  \item{kmer}{the k-mer sequence.}
  \item{position}{the position in the read.}
  \item{count}{the frequency of the k-mer at this position.}
}

\section{Methods}{
  \describe{
    \item{\code{signature(x="SequenceSummary")}}{
      \code{getKmer} is an accessor function that only works if there is
      k-mer data, thus it only works if \code{readSeqFile} was called
      with \code{kmer=TRUE} (and \code{hash.prop} is greater than 0).}
  }
}

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
  ## Load a FASTQ file, with sequence and k-mer hashing on by default.
  s.fastq <- readSeqFile(system.file('extdata', 'test.fastq', package='qrqc'))

  ## plot counts of a subset of k-mers by position
  s.kmers <- getKmer(s.fastq)
  top.kmers <- s.kmers$kmer[order(s.kmers$count, decreasing=TRUE)[1:40]]
  p <- ggplot(subset(s.kmers, kmer \%in\% top.kmers)) + geom_bar(aes(x=position, y=count,
    fill=kmer), stat="identity")
  p
}
\seealso{getGC,getSeqlen,getBase,getBaseProp,getMCQual,kmerKLPlot,kmerEntropyPlot}
\keyword{methods}
\keyword{accessor}
