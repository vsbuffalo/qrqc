\name{getSeqlen-methods}
\docType{methods}
\alias{getSeqlen}
\alias{getSeqlen-methods}
\alias{getSeqlen,SequenceSummary-method}
\title{Get a Data Frame of Sequence Lengths from a \code{SequenceSummary} object}
\description{
  An object that inherits from class \code{SequenceSummary} contains
  sequence length data by position gathered by \code{readSeqFile}. \code{getSeqlen}
  is an accessor function that returns the sequence length data.

  This accessor function is useful if you want to map variables to
  custom \code{ggplot2} aesthetics.
}


\usage{
  getSeqlen(x)
}

\arguments{
  \item{x}{an S4 object that inherits from \code{SequenceSummary} from
    \code{readSeqFile}.}
}


\value{
  \code{getSeqlen} returns a \code{data.frame} with columns:

  \item{length}{the sequence length.}
  \item{count}{the number of reads with this sequence length.}
}

\section{Methods}{
  \describe{
    \item{\code{signature(x = "SequenceSummary")}}{
      \code{getSeqlen} is an accessor function that works on any object read
      in with \code{readSeqFile}; that is, objects that inherit from
      \code{SequenceSummary}.
    }
}}

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
  library(ggplot2)

  ## Load a FASTQ file, with sequence hashing.
  s.trimmed.fastq <- readSeqFile(system.file('extdata', 'test-trimmed.fastq',
    package='qrqc'))

  # A custom plot - a bit contrived, but should show power
  d <- merge(getSeqlen(s.trimmed.fastq), getQual(s.trimmed.fastq),
    by.x="length", by.y="position")
  ggplot(d) + geom_linerange(aes(x=length, ymin=0, ymax=count),
    color="grey") + geom_linerange(aes(x=length, ymin=lower, ymax=upper),
    color="blue") + scale_y_continuous("quality/count") + theme_bw()
}

\seealso{\code{\link{getGC}}, \code{\link{getBase}},
  \code{\link{getBaseProp}}, \code{\link{getQual}},
  \code{\link{getMCQual}}, \code{\link{seqlenPlot}}}

\keyword{methods}
\keyword{accessor}
