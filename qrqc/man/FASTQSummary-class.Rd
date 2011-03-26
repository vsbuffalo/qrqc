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
}
\section{Slots}{
  In addition to the slots inherited from
  \code{\linkS4class{SequenceSummary}},
  \code{\linkS4class{FASTQSummary}} contains:
  
  \describe{ \item{\code{quality}}{a string indicating the type of
    quality (used to convert ASCII characters to quality
    integers). Either "phred", "solexa", or "illumina".}

    \item{\code{qual.freqs}}{a dataframe of quality frequencies by
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
}

\keyword{classes}