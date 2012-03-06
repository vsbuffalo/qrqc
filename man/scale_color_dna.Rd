\name{scale_color_dna}
\alias{scale_color_dna}
\title{Set the color scheme to \code{biovizBase}'s for DNA}
\description{
  This wraps \code{ggplot2}'s \code{scale_color_manual} to use
  \code{biovizBase}'s scheme for DNA (with N).
}

\usage{
scale_color_dna()
}

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
  ## Load a FASTQ file
  s.fastq <- readSeqFile(system.file('extdata', 'test.fastq', package='qrqc'))

  ## Plot qualities with the DNA color scheme
  ggplot(getBase(s.fastq)) + geom_line(aes(x=position, y=frequency,
    color=base)) + scale_color_dna()
}

\seealso{scale_color_iupac,plotBase}
\keyword{graphics}
