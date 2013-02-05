## knitr-method.R -- new report generation using knitr
PDFLATEX.CMD <- "pdflatex -halt-on-error"
TEMPLATES <- list(PDF=system.file('extdata', 'knitr-report-template.Rnw', package='qrqc'))
FIG.HEIGHT <- 5
FIG.WIDTH <- 6

reportName <- function(x) {
  tmp <- sapply(x, function(y) strsplit(basename(y@filename), '.', fixed=TRUE)[[1]][1])
  return(paste(tmp, sep="-"))
}

## setMethod("makeReport", signature(x="FASTQSummary"),
## # Given a FASTQ object, knit a report          
## function(x, output=c("pdf", "html")) {
##   knit(TEMPLATES$PDF)
## })

fileTable <-
function(x) {
  tbl.list <- lapply(x, function(y) {
    n.seqs <- sum(y@seq.lengths)
    seq.lens <- seqLengthRange(y)
    if (y@hashed)
      prop.unique <- length(y@hash)/n.seqs
    list(filename=basename(y@filename), quality=y@quality,
         min.len=seq.lens[1], max.len=seq.lens[2],
         total=n.seqs, prop.unique=round(prop.unique, 3))
  })
  tbl <- do.call(rbind, tbl.list)
  return(xtable(tbl))
}

report <-
# report generates a report using knitr. It accepts a list of objects
# or a single object that inherit from SequenceSummary. In the latter
# case, we make a single list of it; this allows downstream functions
# to work on a consistent object class. Using a generic method is not
# really needed.
function(x, outdir=".", landscape=FALSE) {
  ll <- list()
  if (is(x, "SequenceSummary"))
    ll[[x@filename]] <- x
  else if (is.list(x))
    ll <- x
  else
    stop("report() only handles lists and objects that inherit from SequenceSummary")

  if (landscape) {
    message("using landscape layout")
    FIG.HEIGHT <- 8
    FIG.WIDTH <- 11
  }
  include.hash <- any(sapply(x, function(x) x@hashed))
  include.kmers <- any(sapply(x, function(x) x@kmers.hashed))
  outfile <- file.path(outdir, paste(reportName(x), "tex", sep="."))
  knit(TEMPLATES$PDF)
}

makeHashTable <- function(x, n=16) {
  d <- x@hash[1:n]
  tbl <- as.table(cbind(sequence=names(d), count=d, 'proportion of total'=d/sum(x@seq.lengths)))
  rownames(tbl) <- NULL
  x.tbl <- xtable::xtable(tbl)
  return(x.tbl)
}
