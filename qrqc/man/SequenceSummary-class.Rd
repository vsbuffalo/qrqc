\name{SequenceSummary-class}
\docType{class}
\alias{SequenceSummary-class}
\alias{SequenceSummary}
\title{"SequenceSummary" class representing the summaries of a sequence file}
\description{
  A sequence file read in with \code{readSeqFile} is summarized by a C
  call. Summary statistics on the sequence lengths and base frequencies
  by position (cycle) are stored in this class. If the file is a FASTQ
  file, quality distribution information is also stored. If \code{hash}
  was \code{TRUE} when \code{readSeqFile} was called, the hashed
  sequences are also stored.}
\section{Slots}{
  \code{SequenceSummary} has variable slots, depending on whether the
  sequence file contained quality information, and whetherh the unique
  sequences were hashed.
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
    \item{\code{hash}}{a vector of the frequencies of sequencies (given
      by the name).}
    \item{\code{hashed}}{a logical indicating whether the sequences were
    hased in \code{readSeqFile}.}
  }
}

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
showClass("SequenceSummary")
}

\keyword{classes}