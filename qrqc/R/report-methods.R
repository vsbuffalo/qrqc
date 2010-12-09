## report-methods.R - Methods to create Sweave reports for PDF output of quality information

reportName <- function(obj)
  sprintf("%s-report.html", strsplit(basename(obj@filename), '.', fixed=TRUE)[[1]][1])

setMethod(makeReport, "FASTASummary",
          function(obj) {
            if (!file.exists('images'))
              dir.create('images')
            type <- "FASTA"
            brew(system.file('extdata/report-template.html', package='qrqc'), output=reportName(obj))
          })

setMethod(makeReport, "FASTQSummary",
          function(obj) {
            if (!file.exists('images'))
              dir.create('images')
            type <- "FASTQ"
            html.file <- brew(system.file('extdata/report-template.html', package='qrqc'))
            brew(system.file('extdata/report-template.html', package='qrqc'), output=reportName(obj))
          })

