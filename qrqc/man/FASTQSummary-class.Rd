\name{FASTQSummary-class}
\docType{class}
\alias{FASTQSummary-class}
\alias{FASTQSummary}
\title{"FASTQSummary" class representing the summaries of a FASTQ file}
\description{
  This class contains the same slots as the
  \code{\linkS4class{SequenceSummary}}, as well as additional slots for quality
  information.
}
\section{Slots}{
  \describe{
    \item{\code{base.freqs}}{a dataframe of base frequencies by
      position. Each column is a nucleotide (there is a column for
      position too), and each row contains the count frequencies of bases for that
      position.}
    \item{\code{base.prop}}{a dataframe of base proportions by
      position. This dataframe is stacked, so each row contains a
      single nucleotide, its position, and its proportion for that
      position.}
    \item{\code{seq.lengths}}{a numeric vector of the number of sequences of a
      particular length (the length is the position in the vector).}
    \item{\code{hash}}{a numeric vector of the count frequencies of
      sequences (the sequences are in the name attribute).}
    \item{\code{hashed}}{a logical indicating whether the sequences were
      hashed in \code{readSeqFile}.}
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

\keyword{classes}