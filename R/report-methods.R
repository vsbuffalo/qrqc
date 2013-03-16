## report-methods.R - Methods to create reports for PDF output of quality information
WIDTH <- 8
HEIGHT <- 6.5
DPI <- 72

getReportName <- function(x) {
  .Deprecated("reportName")
  strsplit(basename(x@filename), '.', fixed=TRUE)[[1]][1]
}

makeHashTable <- function(x, n=16) {
  .Deprecated()
  d <- x@hash[1:n]
  tbl <- as.table(cbind(sequence=names(d), count=d, 'proportion of total'=d/sum(x@seq.lengths)))
  rownames(tbl) <- NULL
  x.tbl <- xtable::xtable(tbl)
  return(x.tbl)
}

setMethod(makeReport, "FASTASummary",
          function(x, outputDir=".") {
            .Deprecated("report")
            type <- "FASTA"
            sl.range <- seqLengthRange(x)
            dir <- file.path(makeReportDir(x, outputDir))
            brew::brew(system.file('extdata', 'fasta-report-template.html', package='qrqc'),
                 output=file.path(dir, "report.html"))
            message(sprintf("Report written to directory '%s'.", dir))
          })

setMethod(makeReport, "FASTQSummary",
          function(x, outputDir=".") {
            .Deprecated("report")
            type <- "FASTQ"
            sl.range <- seqLengthRange(x)
            dir <- file.path(makeReportDir(x, outputDir))
            brew::brew(system.file('extdata', 'fastq-report-template.html', package='qrqc'),
                 output=file.path(dir, "report.html"))
            message(sprintf("Report written to directory '%s'.", dir))
          })


