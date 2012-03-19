\name{gcPlot-methods}
\docType{methods}
\alias{gcPlot}
\alias{gcPlot-methods}
\alias{gcPlot,list-method}
\alias{gcPlot,SequenceSummary-method}
\title{Plot GC Content by Position}
\description{
  \code{gcPlot} plots the GC content by position in the read.
}

\usage{
  gcPlot(x, color="red")
}

\arguments{
  \item{x}{an S4 object that inherits from \code{SequenceSummary} from
    \code{readSeqFile}.}
  \item{color}{the color to use for the GC content line.}
}

\section{Methods}{
\describe{
  \item{\code{signature(x = "FASTQSummary")}}{
    \code{gcPlot} will plot the GC content for a single object that
    inherits from \code{SequenceSummary}.
  }
  
  \item{\code{signature(x = "list")}}{
    \code{gcPlot} will plot the GC content for each of the
    \code{SequenceSummary} items in the list and display them in a
    series of panels.
  }
}}

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
  ## Load a FASTQ file, with sequence hashing.
  s.fastq <- readSeqFile(system.file('extdata', 'test.fastq', package='qrqc'))

  ## Plot GC content
  gcPlot(s.fastq)

  ## Plot multiple GC content plots
  s.trimmed.fastq <- readSeqFile(system.file('extdata',
    'test-trimmed.fastq', package='qrqc'))
  gcPlot(list("not trimmed"=s.fastq, "trimmed"=s.trimmed.fastq))

  ## Graphical features can be added
  gcPlot(s.trimmed.fastq) + geom_hline(yintercept=0.5, color="purple")
}

\seealso{\code{\link{getBase}}, \code{\link{getBaseProp}}}

\keyword{methods}
\keyword{graphics}

