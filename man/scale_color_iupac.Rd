\name{scale_color_iupac}
\alias{scale_color_iupac}
\title{Set the color scheme to \code{biovizBase}'s for IUPAC codes}
\description{
  This wraps \code{ggplot2}'s \code{scale_color_manual} to use
  \code{biovizBase}'s scheme IUPAC nucleotides codes.
}

\usage{
scale_color_iupac()
}

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
  ## Load a FASTQ file
  s.fastq <- readSeqFile(system.file('extdata', 'test.fastq', package='qrqc'))

  ## Plot qualities with the DNA color scheme
  ggplot(getBase(s.fastq)) + geom_line(aes(x=position, y=frequency,
    color=base)) + scale_color_iupac()
}

\seealso{scale_color_dna,plotBase}
\keyword{graphics}
