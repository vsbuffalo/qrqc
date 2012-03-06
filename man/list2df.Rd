\name{list2df}
\alias{list2df}
\title{Apply a function to items in list and combine into data frame}
\description{
  \code{list2df} is a helper function that takes a named list and
  applies a function to each element, and combines the resulting data
  frames into a single data frame. The output data frame will have an
  additional column named \code{sample} indicating which element the
  data came from.
}

\usage{
  list2df(x, fun)
}

\arguments{
  \item{x}{a named list of objects.}
  \item{fun}{a function that takes in the elements of \code{x} and
  outputs a data frame.}
}

\value{
  A data frame made by applying \code{fun} to each element of the list
  \code{x}. An additional column named \code{sample} will indicate which
  element the data came from.
}

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
  ## Get some sequence files
  sq.files = list.files(system.file('extdata', package='qrqc'),
    pattern="test.*fastq", full.names=TRUE)
  names(sq.files) <- gsub("(.*)\\\\.fastq", "\\1", basename(sq.files))
  sq <- lapply(sq.files, readSeqFile)


  ## Take the FASTQSummary objects, extract quality data from each of
  ## the, and combine.
  d <- list2df(sq, getQual)

  ## Look at difference in average quality
  aggregate(d$mean, list(sample=d$sample), mean)

  ## Look at difference in variance - this is where we really see a
  ## change.
  aggregate(d$mean, list(sample=d$sample), var)
}