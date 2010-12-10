## report-methods.R - Methods to create Sweave reports for PDF output of quality information
WIDTH <- 900
HEIGHT <- 500

reportName <- function(obj)
  sprintf("%s-report.html", strsplit(basename(obj@filename), '.', fixed=TRUE)[[1]][1])

seqLengthRange <- function(obj) {
  s <- which(obj@seq.lengths > 0)
  return(c(min(s), max(s)))
}

makeHashTable <- function(obj, n=16) {
  d <- obj@hash[1:n]
  tbl <- as.table(cbind(sequence=names(d), count=d, 'proportion of total'=d/sum(obj@seq.lengths)))
  rownames(tbl) <- NULL
  x.tbl <- xtable(tbl)
  return(x.tbl)
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

