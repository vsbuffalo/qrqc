\name{SequenceSummary-class}
\docType{class}
\alias{SequenceSummary,SequenceSummary-class}
\title{"SequenceSummary" class representing the summaries of a sequence file}
\description{
  A sequence file read in with \code{readSeqFile} is summarized by a C
  call. This is a base class with slots common to both
  \code{\linkS4class{FASTQSummary}} and  \code{\linkS4class{FASTASummary}}. This is
  not usually instantiated directly.
}
\section{Slots}{
  \describe{
    \item{\code{filename}}{the filename processed by \code{readSeqFile}.}
   \item{\code{base.freqs}}{a dataframe of base frequencies by
      position. Each column is a nucleotide (there is a column for
      position too), and each row contains the count frequencies of bases for that
      position.}
    \item{\code{base.props}}{a dataframe of base proportions by
      position. This dataframe is stacked, so each row contains a
      single nucleotide, its position, and its proportion for that
      position.}
    \item{\code{seq.lengths}}{a numeric vector of the number of sequences of a
      particular length (the length is the position in the vector).}
    \item{\code{hash}}{a numeric vector of the count frequencies of
      sequences (the sequences are in the name attribute).}
    \item{\code{hashed}}{a logical indicating whether the sequences were
      hashed in \code{readSeqFile}.}
  }
}

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
showClass("SequenceSummary")
}


\seealso{

  \code{\link[=FASTQSummary-class]{FASTQSummary}} and
  \code{\link[=FASTASummary-class]{FASTASummary}} are the classes that
  inherit from \code{SequenceSummary}.

  \code{\link{readSeqFile}} is the function that takes a FASTQ or FASTA
  file and returns a \code{FASTQSummary} object or \code{FASTASummary}
  object.

}


\keyword{classes}