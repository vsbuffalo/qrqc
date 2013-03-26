\name{report-methods}
\docType{methods}
\alias{report}
\alias{report-methods}
\title{Make a TeX or HTML report from a single FASTQSummary or
  FASTASummary object, or a list of these objects}
\description{
  \code{report} takes single
  \code{\link[=FASTQSummary-class]{FASTQSummary}} or
  \code{\link[=FASTASummary-class]{FASTASummary}} object, or a list of
  these objects and creates a TeX or HTML report, and writes it to a file within a
  directory.
}
\usage{
  report(x, outdir=datedir(), type=c("tex", "html"), kmers=TRUE,
         hash=FALSE, landscape=TRUE, template=NULL, tar=FALSE,
         compression=c("none", "gzip", "bzip2"), quiet=FALSE)
}
\arguments{
  \item{x}{an object that is either a single \code{FASTQSummary} or
    \code{FASTASummary} object, or a list of these objects.}
  \item{outdir}{a character argument to to indicate the report
    output directory. By default, this runs the function
    \code{datedir()} which creates a directory in the format
    "qrqc-report-2013-03-07T212756", where the date is given in ISO 8601
    time.}
  \item{type}{a string indicating the output type: either
    "tex" or "html".}
  \item{kmers}{a logical value indicating whether to include the plot
    create by \code{kmerKLPlot}. Only available if \code{readSeqFile}
    was called with \code{kmer=TRUE}.}
  \item{hash}{a logical value indicating whether to include the table of
    common hashed sequences. Only available if \code{readSeqFile} was
    called with \code{hash=TRUE}.}
  \item{landscape}{a logical value indicating whether to output the TeX
    output in landscape orientation. Only applicable to \code{type="tex"}.}
  \item{template}{an optional alternate knitr Rnw or Rhtml file. If
    \code{NULL}, this will use the default template file. This can also
    be changed with \code{qrqc_options}.}
  \item{tar}{a logical value indicating whether to tar the report
    directory. Optionally, can be compressed via the \code{compression}
    argument.}
  \item{compression}{a string indicating whether to use compression:
    either "none", "gzip", or "bzip2".}
}

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
\dontrun{
  # Not run to avoid cluttering the file system with these files
  ## Load a few FASTQ file
  s.fastq <- readSeqFile(system.file('extdata', 'test.fastq', package='qrqc'))
  s.trimmed.fastq <- readSeqFile("qrqc/inst/extdata/test-trimmed.fastq", kmer=FALSE)

  ## make a TeX report, and tar and compress in
  Report(list(untrimmed=s.fastq), trimmed=s.trimmed.fastq, type="tex",
         tar=TRUE, compression="gzip")
}
}

\keyword{file}