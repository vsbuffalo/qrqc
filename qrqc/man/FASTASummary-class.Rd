\name{FASTASummary-class}
\Rdversion{1.1}
\docType{class}
\alias{FASTASummary-class}
\alias{makeReport,FASTASummary-method}
\title{\code{FASTASummary} class representing the summaries of a FASTA file}
\description{
  This class contains the same slots as the
  \code{\linkS4class{SequenceSummary}}, but it is used to indicate the
  data originated from a FASTA file.

  Note that many accessor functions transform data in the slots into data
  frames. The data in the slots is mostly untransformed and less easy to
  work with directly, so using the accessor functions is recommended.

}
\section{Slots}{
  \code{\linkS4class{FASTASummary}} has the slots inherited from 
  \code{\linkS4class{SequenceSummary}}.
}

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
showClass("FASTASummary")
}

\seealso{

  \code{\link[=FASTQSummary-class]{FASTQSummary}} is the
  counterpart of this class for FASTQ data.

  \code{\link{readSeqFile}} is the function that takes a FASTA file and
  returns a \code{FASTASummary} object.
  
  \code{\link{basePlot}} is a function that plots the distribution of
  bases over sequence length for a particular \code{FASTASummary}
  object. \code{\link{gcPlot}} combines and plots the GC proportion.
  
  \code{\link{seqlenPlot}} is a function that plots a histogram of
  sequence lengths for a particular \code{FASTASummary} object.

  \code{\link{kmerKLPlot}} is a function that uses Kullback-Leibler
  divergence to make a plot that can aid in finding possible
  contamination (if \code{readSeqFile} had \code{kmer=TRUE}).

  \code{\link{kmerEntropyPlot}} is a function that plots the Shannon
  entropy of k-mers per position.

  There are acccessor functions \code{\link{getQual}},
  \code{\link{getBase}}, \code{\link{getBaseProp}},
  \code{\link{getSeqlen}}, \code{\link{getKmer}}, \code{\link{getGC}}
  for transforming the raw data in the object's slot (direct from the C
  call) to more usable data frames.  

}

\keyword{classes}