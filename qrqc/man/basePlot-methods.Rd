\name{basePlot-methods}
\docType{methods}
\alias{basePlot-methods}
\alias{basePlot,list-method}
\alias{basePlot,SequenceSummary-method}
\title{Plot Base Frequency or Proportion by Position}
\description{
  \code{basePlot} plots the frequency or proportion of bases by
  position in the read.

  Specific bases (such as "N") can be plot alone with this function too.
}

\usage{
  basePlot(x, geom=c("line", "bar", "dodge"),
    type=c("frequency", "proportion"), bases=DNA_BASES_N,
    colorvalues=getBioColor("DNA_BASES_N"))
}

\arguments{
  \item{x}{an S4 object that inherits from \code{SequenceSummary} from
    \code{readSeqFile}.}
  \item{geom}{Either "line", "bar", or "dodge" indicating
    the geom to use when plotting the bases. "line" will plot base
    proportion of frequency with lines. "bar" and "dodge" will use bars;
    "bar" defaults to filling the bars with different colors to distinguish
    bases, "dodge" lays the bars side by side.}
  \item{type}{Either "frequency" or "proportion" indicating whether to
    use count data or the proportion per base.}
  \item{bases}{a character vector indicating which bases to include. By
    default, all bases in DNA_BASES_N. Another good option would be
    IUPAC_CODE_MAP, which is included in the \code{Biostrings} package.}
  \item{colorvalues}{a character vectors of colors to use; the names of
    the elements must map to the bases.}
}

\section{Methods}{
\describe{
  \item{\code{signature(x = "FASTQSummary")}}{
    \code{basePlot} will plot the base frequencies or proportions for a
    single object that inherits from \code{SequenceSummary}.   
  }
  
  \item{\code{signature(x = "list")}}{
    \code{basePlot} will plot the base frequencies or proportions for each of the
    \code{SequenceSummary} items in the list and display them in a series of
    panels.}
}}

\author{Vince Buffalo <vsbuffalo@ucdavis.edu>}
\examples{
  ## Load a FASTQ file, with sequence hashing.
  s.fastq <- readSeqFile(system.file('extdata', 'test.fastq', package='qrqc'))

  ## Plot bases
  basePlot(s.fastq)

  ## Plot bases with filled bars
  basePlot(s.fastq, geom="bar")

  ## Plot bases with dodged bars
  basePlot(s.fastq, geom="dodge")

  ## Plot bases with dodged bars
  basePlot(s.fastq, geom="bar", bases=c("G", "T"))

  ## Plot multiple base plots
  s.trimmed.fastq <- readSeqFile(system.file('extdata',
    'test-trimmed.fastq', package='qrqc'))
  basePlot(list("not trimmed"=s.fastq, "trimmed"=s.trimmed.fastq))


  ## Graphical features can be added
  basePlot(s.trimmed.fastq, type="proportion") +
    geom_hline(yintercept=0.25, color="purple")
}

\seealso{getBase}
\seealso{getBaseProp}
\keyword{methods}
\keyword{graphics}

