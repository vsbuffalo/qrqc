\name{FASTASummary-class}
\docType{class}
\alias{FASTASummary,FASTASummary-class}
\title{"FASTASummary" class representing the summaries of a FASTA file}
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
}

\keyword{classes}