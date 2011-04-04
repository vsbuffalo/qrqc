\name{readSeqFile}
\alias{readSeqFile}
\title{Read and Summarize a Sequence (FASTA or FASTQ) File}
\description{
  \code{readSeqFile} reads a FASTQ or FASTA file, summarizing the
  nucleotide distribution across position (cycles) and the sequence
  length distributions. If \code{type} is `fastq', the distribution
  of qualities across position will also be recorded. If \code{hash} is
  \code{TRUE}, the unique sequences will be hashed with counts of their
  frequency.
}
\usage{
  readSeqFile(filename, type='fastq', max.length=1000,
              quality='illumina', hash=TRUE, verbose=FALSE)
}
\arguments{
  \item{filename}{the name of the file which the sequences are to be
    read from.}
  \item{type}{either `fastq' or `fasta', representing the type of the
    file. FASTQ files will have the quality distribution by position
    summarized.}
  \item{max.length}{the largest sequence length likely to be
    encountered. For efficiency, a matrix larger than the largest
    sequence is allocated to *this* size in C, populated, and then
    trimmed in R. Specifying a value too small will lead to an error and
    the function will need to be re-run.}
  \item{quality}{either `illumina', `phred', or `solexa', this
    determines the quality offsets and range. See the values of
    QUALITY.CONSTANTS for more information.}
  \item{hash}{a logical value indicating whether to hash sequences}
  \item{verbose}{a logical value indicating whether be verbose (in the C
    backend).}
}
\value{
  An S4 object of \code{\linkS4class{FASTQSummary}} or
  \code{\linkS4class{FASTASummary}} containing the summary statistics.
}
\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
  ## Load a FASTQ file, with sequence hashing.
  s.fastq <- readSeqFile(system.file('extdata', 'test.fastq', package='qrqc'))

  ## Load a FASTA file, without sequence hashing.
  s.fasta <- readSeqFile(system.file('extdata', 'test.fasta', package='qrqc'),
                         type='fasta', hash=FALSE)

}

\seealso{ \code{\link[=FASTQSummary-class]{FASTQSummary}} and
  \code{\link[=FASTASummary-class]{FASTASummary}} are the classes used
  by the objects returned by \code{readSeqFile}.}

\keyword{file}