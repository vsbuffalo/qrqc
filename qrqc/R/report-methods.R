## report-methods.R - Methods to create reports for PDF output of quality information
WIDTH <- 8
HEIGHT <- 6.5
DPI <- 72

getReportName <- function(x)
  strsplit(basename(x@filename), '.', fixed=TRUE)[[1]][1]

seqLengthRange <- function(x) {
  s <- which(x@seq.lengths > 0)
  return(c(min(s), max(s)))
}

makeHashTable <- function(x, n=16) {
  d <- x@hash[1:n]
  tbl <- as.table(cbind(sequence=names(d), count=d, 'proportion of total'=d/sum(x@seq.lengths)))
  rownames(tbl) <- NULL
  x.tbl <- xtable::xtable(tbl)
  return(x.tbl)
}

makeReportDir <- function(x, outputDir=".") {
  current.reports <- dir(outputDir, pattern=".*-report")
  if (length(current.reports)) {
    # we need to find the current report number and create a directory
    # with the digit incremented.

    last <- 0
    if (length(grep(".*-report$", current.reports))) {
      # There's a report is first, but will not match -x where x is an
      # integer.
      last <- 1
    }
    
    tmp <- gsub(".*-report-(\\d+)", '\\1', current.reports)
    last <- suppressWarnings(max(na.exclude(c(as.numeric(tmp), last))))
    dirpath <- file.path(outputDir, sprintf("%s-report-%s", getReportName(x), last+1))
  } else {
    dirpath <- file.path(outputDir, sprintf("%s-report", getReportName(x)))
  }  
  if (!(dir.create(dirpath) && dir.create(file.path(dirpath, "images"))))
    stop("Could not create report directory; perhaps permissions do not allow this.")
  return(dirpath)
}


setMethod(makeReport, "FASTASummary",
          function(x, outputDir=".") {
            type <- "FASTA"
            sl.range <- seqLengthRange(x)
            dir <- file.path(makeReportDir(x, outputDir))
            brew::brew(system.file('extdata', 'fasta-report-template.html', package='qrqc'),
                 output=file.path(dir, "report.html"))
            message(sprintf("Report written to directory '%s'.", dir))
          })

setMethod(makeReport, "FASTQSummary",
          function(x, outputDir=".") {
            type <- "FASTQ"
            sl.range <- seqLengthRange(x)
            dir <- file.path(makeReportDir(x, outputDir))
            brew::brew(system.file('extdata', 'fastq-report-template.html', package='qrqc'),
                 output=file.path(dir, "report.html"))
            message(sprintf("Report written to directory '%s'.", dir))
          })


