\name{readSeqFile}
\alias{readSeqFile}
\title{Read and Summarize a Sequence (FASTA or FASTQ) File}
\description{
  \code{readSeqFile} reads a FASTQ or FASTA file, summarizing the
  nucleotide distribution across position (cycles) and the sequence
  length distributions. If \code{type} is `fastq', the distribution
  of qualities across position will also be recorded. If \code{hash} is
  \code{TRUE}, the unique sequences will be hashed with counts of their
  frequency. By default, only 10\% of the reads will be hashed; this
  proportion can be controlled with \code{hash.prop}. If
  \code{kmer=TRUE}, k-mers of length \code{k} will be hashed by
  position, also with the sampling proportion controlled by
  \code{hash.prop}.
  
  
}
\usage{
  readSeqFile(filename, type=c("fastq", "fasta"), max.length=1000,
              quality=c("sanger", "solexa", "illumina"), hash=TRUE,
              hash.prop=0.1, kmer=TRUE, k=6L, verbose=FALSE)
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
  \item{quality}{either `illumina', `sanger', or `solexa', this
    determines the quality offsets and range. See the values of
    QUALITY.CONSTANTS for more information.}
  \item{hash}{a logical value indicating whether to hash sequences}
  \item{hash.prop}{a numeric value in (0, 1] that functions as the proportion
    of reads to hash.}
  \item{kmer}{a logical value indicating whether to hash k-mers by
    position.}
  \item{k}{an integer value indicating the k-mer size.}
  \item{verbose}{a logical value indicating whether be verbose (in the C
    backend).}
}
\value{
  An S4 object of \code{\linkS4class{FASTQSummary}} or
  \code{\linkS4class{FASTASummary}} containing the summary statistics.
}

\note{
  Identifying the correct quality can be difficult. \code{readSeqFile}
  will error out if it a base quality outside of the range of a known
  quality type, but it is possible one could have reads with a different
  quality type that won't fall outside of the another type.

  Here is a bit more about quality:

  \item{phred}{PHRED quality scores (e.g. from Roche 454). ASCII with no offset,
    range: [4, 60]. This has been removed as an option since sequence reads with this
    type are very, very uncommon.}
  \item{sanger}{Sanger are PHRED ASCII qualities with an offset of 33, range:
    [0, 93]. From NCBI SRA, or Illumina pipeline 1.8+.}
  \item{solexa}{Solexa (also very early Illumina - pipeline < 1.3). ASCII offset
    of 64, range: [-5, 62]. Uses a different quality-to-probabilities
    conversion than other schemes.}
  \item{illumina}{Illumina output from pipeline versions between 1.3 and
    1.7. ASCII offset of 64, range: [0, 62].}
}

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
  ## Load a FASTQ file, with sequence hashing.
  s.fastq <- readSeqFile(system.file('extdata', 'test.fastq', package='qrqc'))

  ## Load a FASTA file, without sequence hashing.
  s.fasta <- readSeqFile(system.file('extdata', 'test.fasta', package='qrqc'),
                         type='fasta', hash=FALSE)

}

\seealso{

  \code{\link[=FASTQSummary-class]{FASTQSummary}} and
  \code{\link[=FASTASummary-class]{FASTASummary}} are the classes of the
  objects returned by \code{readSeqFile}.

  \code{\link{plotBases}} is a function that plots the distribution of
  bases over sequence length for a particular \code{FASTASummary} or
  \code{FASTQSummary} object. \code{\link{plotGC}} combines and plots
  the GC proportion.
  
    
  \code{\link{plotQuals}} is a function that plots the distribution of
  qualities over sequence length for a particular \code{FASTASummary}
  or \code{FASTQSummary} object.

  \code{\link{plotSeqLengths}} is a function that plots a histogram of
  sequence lengths for a particular \code{FASTASummary} or
  \code{FASTQSummary} object.

}



\keyword{file}