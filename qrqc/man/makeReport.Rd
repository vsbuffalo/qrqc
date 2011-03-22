\name{makeReport}
\alias{makeReport,FASTASummary-method}
\title{Make an HTML report from a FASTASummary of FASTQSummary object}
\description{
  \code{makeReport} takes a \code{FASTQSummary} or \code{FASTASummary}
  object, creates an HTML report, and writes it to a file.
}
\usage{
  makeReport(obj, filename=NULL)
}
\arguments{
  \item{obj}{an object that is either \code{FASTQSummary} or \code{FASTASummary}.}
  \item{filename}{optional string that is the filename that the report
    will be written to. By default, this is \code{NULL}, indicating the
    file will be saved by default under _filename_-report.html (where
    _filename_ is from the object's \code{filename} slot.)}
}
\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
  ## Load a FASTQ file, with sequence hashing.
  s.fastq <- readSeqFile(system.file('extdata', 'test.fastq', package='qrqc'))

  ## Make and save a report
  makeReport(s.fastq)
}

\keyword{file}