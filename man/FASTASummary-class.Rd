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
  
  \code{\link{plotBases}} is a function that plots the distribution of
  bases over sequence length for a particular \code{FASTASummary}
  object. \code{\link{plotGC}} combines and plots the GC proportion.
  
  \code{\link{plotSeqLengths}} is a function that plots a histogram of
  sequence lengths for a particular \code{FASTASummary} object.

}

\keyword{classes}