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
    \item{\code{base.freqs}}{a dataframe of base frequencies by position
      (cycle). Each column is a nucleotide (and there's a column for
      position), and each row contains the frequencies of bases for that
      position.}
    \item{\code{base.prop}}{a dataframe of base proportions by position
      (cycle). This dataframe is unstacked, so each row contains a
      single nucleotide, its position, and its proportion for that
      position.}
    \item{\code{qual.freqs}}{a dataframe of quality frequencies by
      position, if the file was a FASTQ file.}
    \item{\code{mean.qual}}{a numeric that is the mean quality across all positions
      (cycles), weighted by the number of reads that extended to that
      position.}
    \item{\code{seq.lengths}}{a vector of the number of sequences of a
      particular length (the position in the vector).}
    \item{\code{quality}}{a string representing the quality scheme.}
    \item{\code{type}}{a string representing the type of file input;
      either `FASTQ' or `FASTA'.}
    \item{\code{hash}}{a vector of the frequencies of sequences.}
    \item{\code{hashed}}{a logical indicating whether the sequences were
      hased in \code{readSeqFile}.}
  }
}

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
showClass("FASTQSummary")
}

\keyword{classes}