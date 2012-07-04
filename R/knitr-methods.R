## knitr-method.R -- new report generation using knitr

getReportName <- function(x) {
  strsplit(basename(x@filename), '.', fixed=TRUE)[[1]][1]
}

## setMethod("makeReport", signature(x="FASTQSummary"),
## # Given a FASTQ object, knit a report          
## function(x, output=c("pdf", "html")) {
##   stitch()
## })


