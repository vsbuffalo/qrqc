\name{makeReport-methods}
\docType{methods}
\alias{makeReport}
\alias{makeReport-methods}
\title{Make an HTML report from a FASTASummary of FASTQSummary object}
\description{
  \code{makeReport} takes a \code{FASTQSummary} or \code{FASTASummary}
  object, creates an HTML report, and writes it to a file within a
  directory. The directory naming is incremental so past reports will
  not be overwritten.
}
\usage{
  makeReport(x, outputDir=".")
}
\arguments{
  \item{x}{an object that is either \code{FASTQSummary} or \code{FASTASummary}.}
  \item{outputDir}{an optional character argument to to indicate the report
    output directory. By default, the current directory.}
}
\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
  ## Load a FASTQ file
  s.fastq <- readSeqFile(system.file('extdata', 'test.fastq', package='qrqc'))

  ## Make and save a report
  makeReport(s.fastq)
}

\keyword{file}