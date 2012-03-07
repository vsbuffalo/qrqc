\name{getBaseProp-methods}
\docType{methods}
\alias{getBaseProp}
\alias{getBaseProp-methods}
\alias{getBaseProp,SequenceSummary-method}
\title{Get a Data Frame of Base Proportion Data from a \code{SequenceSummary} object}
\description{
  An object that inherits from class \code{SequenceSummary} contains
  base frequency data by position gathered by \code{readSeqFile}. \code{getBaseProp}
  is an accessor function that reshapes the base frequency data by position
  into a data frame and calculates the proportions of each base per position.

  This accessor function is useful if you want to map variables to
  custom \code{ggplot2} aesthetics. Base frequency be accessed
  with \code{getBase}.
}


\usage{
  getBaseProp(x, drop=TRUE)
}

\arguments{
  \item{x}{an S4 object that inherits from \code{SequenceSummary} from
    \code{readSeqFile}.}
  \item{drop}{a logical value indicating whether to drop bases that
  don't have any counts.}
}


\value{
  \code{getBaseProp} returns a \code{data.frame} with columns:

  \item{position}{the position in the read.}
  \item{base}{the base.}
  \item{proportion}{the proportion of a base found per position in the read.}
}

\section{Methods}{
  \describe{
    \item{\code{signature(x = "SequenceSummary")}}{
      \code{getBaseProp} is an accessor function that works on any object read
      in with \code{readSeqFile}; that is, objects that inherit from
      \code{SequenceSummary}.
    }
}}

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
  ## Load a FASTQ file, with sequence hashing.
  s.fastq <- readSeqFile(system.file('extdata', 'test.fastq',
    package='qrqc'))

  # A custom base plot
  ggplot(getBaseProp(s.fastq)) + geom_line(aes(x=position, y=proportion,
    color=base)) + facet_grid(. ~ base) + scale_color_dna()
}

\seealso{getGC,getSeqlen,getBase,getQual,getMCQual,basePlot}
\keyword{methods}
\keyword{accessor}
