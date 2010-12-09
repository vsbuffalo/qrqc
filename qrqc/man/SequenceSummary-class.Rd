\name{SequenceSummary-class}
\docType{class}
\alias{SequenceSummary-class}
\alias{SequenceSummary}
\title{"SequenceSummary" class representing the summaries of a sequence file}
\description{
  A sequence file read in with \code{readSeqFile} is summarized by a C
  call. This is a base class with slots common to both
  \code{\linkS4class{FASTQSummary}} and  \code{\linkS4class{FASTASummary}}. This is
  not usually instantiated directly.
}
\section{Slots}{
  \describe{
    \item{\code{base.freqs}}{a dataframe of base frequencies by position
      (cycle). Each column is a nucleotide (and there's a column for
      position), and each row contains the frequencies of bases for that
      position.}
    \item{\code{base.prop}}{a dataframe of base proportions by position
      (cycle). This dataframe is unstacked, so each row contains a
      single nucleotide, its position, and its proportion for that
      position.}
    \item{\code{seq.lengths}}{a vector of the number of sequences of a
      particular length (the position in the vector).}
    \item{\code{quality}}{a string representing the quality scheme.}
    \item{\code{hash}}{a vector of the frequencies of sequences.}
    \item{\code{hashed}}{a logical indicating whether the sequences were
    hased in \code{readSeqFile}.}
  }
}

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
showClass("SequenceSummary")
}

\keyword{classes}