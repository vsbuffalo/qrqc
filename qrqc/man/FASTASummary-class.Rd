\name{FASTASummary-class}
\docType{class}
\alias{FASTASummary-class}
\alias{FASTASummary}
\title{"FASTASummary" class representing the summaries of a FASTA file}
\description{
  This class contains the same slots as the
  \code{\linkS4class{SequenceSummary}}, but is used to indicate a sequence file
  is a FASTA-formatted file.
}
\section{Slots}{
  \describe{
    \item{\code{base.freqs}}{a dataframe of base frequencies by
      position. 
      Each column is a nucleotide (and there's a column for
      position), and each row contains the count frequencies of bases for that
      position.}
    \item{\code{base.prop}}{a dataframe of base proportions by
      position. 
      This dataframe is unstacked, so each row contains a
      single nucleotide, its position, and its proportion for that
      position.}
    \item{\code{seq.lengths}}{a vector of the number of sequences of a
      particular length (the position in the vector).}
    \item{\code{hash}}{a vector of the frequencies of sequences.}
    \item{\code{hashed}}{a logical indicating whether the sequences were
      hased in \code{readSeqFile}.}
}
}

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
showClass("FASTASummary")
}

\keyword{classes}