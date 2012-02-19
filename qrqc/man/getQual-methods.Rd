\name{getQual-methods}
\docType{methods}
\alias{getQual-methods}
\alias{getQual,FASTQSummary-method}
\title{Get a Data Frame of Quality Data from a \code{FASTQSummary} object}
\description{
  An object of class \code{FASTQSummary} contains quality data (binned
  by \code{readSeqFile}). \code{getQual} is an accessor function that
  reshapes the data into a data frame.

  This accessor function is useful if you want to map variables to
  custom \code{ggplot2} aesthetics.
}


\usage{
  getQual(x)
}

\arguments{
  \item{x}{an S4 object of class \code{FASTQSummary} from
    \code{readSeqFile}.}
}

\value{
  \code{qualPlot} returns a \code{data.frame} with columns:

  \item{position}{the position in the read.}
  \item{ymin}{the minimum quality found per a position in the read.}
  \item{alt.lower}{the 10\% quantile found per a position in the read.}
  \item{lower}{the 25\% quartile found per a position in the read.}
  \item{middle}{the median found per a position in the read.}
  \item{upper}{the 75\% quartile found per a position in the read.}
  \item{alt.upper}{the 90\% quantile found per a position in the read.}
  \item{ymax}{the maximum quality found per a position in the read.}
  \item{mean}{the mean quality (calculated from the binned data by using a
  weighted mean function) per the position in the read.}
}

\section{Methods}{
\describe{
\item{\code{signature(x = "FASTQSummary")}}{
  \code{getQual} is an accessor function that only works if there is
  quality data, thus it only works with objects of class
  \code{FASTQSummary}.   
}
}}

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
  ## Load a FASTQ file, with sequence hashing.
  s.fastq <- readSeqFile(system.file('extdata', 'test.fastq',
    package='qrqc'))

  ## Mean quality by position
  p <- ggplot(getQual(s.fastq)) + geom_line(aes(x=position, y=mean),
    color="blue")
  p
  p <- p + scale_y_continuous(limits=c(0, 42))
  p
}
\seealso{getQual}
\seealso{list2df}
\keyword{methods}
\keyword{accessor}
