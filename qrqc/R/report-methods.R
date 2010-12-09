## report-methods.R - Methods to create Sweave reports for PDF output of quality information
WIDTH <- 900
HEIGHT <- 500

reportName <- function(obj)
  sprintf("%s-report.html", strsplit(basename(obj@filename), '.', fixed=TRUE)[[1]][1])

seqLengthRange <- function(obj) {
  s <- which(obj@seq.lengths > 0)
  return(c(min(s), max(s)))
}
  
setMethod(makeReport, "FASTASummary",
          function(obj) {
            if (!file.exists('images'))
              dir.create('images')
            type <- "FASTA"
            sl.range <- seqLengthRange(obj)
            brew(system.file('extdata/fasta-report-template.html', package='qrqc'), output=reportName(obj))
          })

setMethod(makeReport, "FASTQSummary",
          function(obj) {
            if (!file.exists('images'))
              dir.create('images')
            type <- "FASTQ"
            sl.range <- seqLengthRange(obj)
            brew(system.file('extdata/fastq-report-template.html', package='qrqc'), output=reportName(obj))
          })

