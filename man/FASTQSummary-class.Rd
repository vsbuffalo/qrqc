\name{FASTQSummary-class}
\Rdversion{1.1}
\docType{class}
\alias{FASTQSummary-class}
\alias{makeReport,FASTQSummary-method}
\alias{plotQuals,FASTQSummary-method}
\title{\code{FASTQSummary} class representing the summaries of a FASTQ file}
\description{
  This class contains the same slots as the
  \code{\linkS4class{SequenceSummary}}, as well as additional slots for quality
  information.

  Note that many accessor functions transform data in the slots into data
  frames. The data in the slots is mostly untransformed and less easy to
  work with directly, so using the accessor functions is recommended.

}
\section{Slots}{
  In addition to the slots inherited from
  \code{\linkS4class{SequenceSummary}},
  \code{\linkS4class{FASTQSummary}} contains:
  
  \describe{ \item{\code{quality}}{a string indicating the type of
    quality (used to convert ASCII characters to quality
    integers). Either "phred", "solexa", or "illumina".}

    \item{\code{qual.freqs}}{a data frame of quality frequencies by
      position, if the file was a FASTQ file.}
    \item{\code{mean.qual}}{a numeric that is the mean quality across
      all positions,
      weighted by the number of reads that extended to that
      position.}
  }
}

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
showClass("FASTQSummary")
}


\seealso{

  \code{\link[=FASTASummary-class]{FASTASummary}} is the
  counterpart of this class for FASTA data.

  \code{\link{readSeqFile}} is the function that takes a FASTQ file and
  returns a \code{FASTQSummary} object.

  \code{\link{basePlot}} is a function that plots the distribution of
  bases over sequence length for a particular \code{FASTQSummary}
  object. \code{\link{gcPlot}} combines and plots the GC proportion.
  
  \code{\link{qualPlot}} is a function that plots the distribution of
  qualities over sequence length for a particular \code{FASTQSummary}
  object.

  \code{\link{seqlenPlot}} is a function that plots a histogram of
  sequence lengths for a particular \code{FASTQSummary} object.

  
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