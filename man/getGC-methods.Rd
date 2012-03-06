\name{getGC-methods}
\docType{methods}
\alias{getGC}
\alias{getGC-methods}
\alias{getGC,SequenceSummary-method}
\title{Get a Data Frame of GC Content from a \code{SequenceSummary} object}
\description{
  An object that inherits from class \code{SequenceSummary} contains
  base frequency data by position gathered by \code{readSeqFile}. \code{getGC}
  is an accessor function that reshapes the base frequency data into a
  data frame and returns the GC content by position.

  This accessor function is useful if you want to map variables to
  custom \code{ggplot2} aesthetics. Frequencies or proportions of all
  bases (not just GC) can be accessed with \code{getBase} and
  \code{getBaseProp} respectively.
}

\usage{
  getGC(x)
}

\arguments{
  \item{x}{an S4 object that inherits from \code{SequenceSummary} from
    \code{readSeqFile}.}
}


\value{
  \code{getGC} returns a \code{data.frame} with columns:

  \item{position}{the position in the read.}
  \item{gc}{GC content per position in the read.}
}

\section{Methods}{
  \describe{
    \item{\code{signature(x = "SequenceSummary")}}{
      \code{getGC} is an accessor function that works on any object read
      in with \code{readSeqFile}; that is, objects that inherit from
      \code{SequenceSummary}.
    }
}}

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
  ## Load a FASTQ file, with sequence hashing.
  s.fastq <- readSeqFile(system.file('extdata', 'test.fastq',
    package='qrqc'))

  # A custom GC plot
  d <- merge(getQual(s.fastq), getGC(s.fastq), by.x="position", by.y="position")
  p <- ggplot(d) + geom_linerange(aes(x=position, ymin=lower,
    ymax=upper, color=gc)) + scale_color_gradient(low="red",
    high="blue") + scale_y_continuous("GC content")
  p
}

\seealso{getSeqlen,getBase,getBaseProp,getQual,getMCQual,cPlot}
\keyword{methods}
\keyword{accessor}
