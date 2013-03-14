## knitr-method.R -- new report generation using knitr
TEMPLATES <- list(PDF=system.file('extdata', 'knitr-report-template.Rnw', package='qrqc'))
FIG.HEIGHT <- 6
FIG.WIDTH <- 8

reportName <- function(x) {
  if (is.list(x)) {
    tmp <- sapply(x, function(y) strsplit(basename(y@filename), '.', fixed=TRUE)[[1]][1])
    out <- paste(tmp, sep="-")
  } else {
    out <- basename(y@filename)
  }
  return(out)
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

fileTable <-
function(x) {
  tbl.list <- lapply(x, function(y) {
    n.seqs <- sum(y@seq.lengths)
    seq.lens <- seqLengthRange(y)
    if (y@hashed)
      prop.unique <- length(y@hash)/n.seqs
    list(Filename=basename(y@filename), Quality=y@quality,
         "Min Length"=seq.lens[1], "Max Length"=seq.lens[2],
         Total=n.seqs, "Proportion Unique"=round(prop.unique, 3))
  })
  tbl <- do.call(rbind, tbl.list)
  tbl <- as.data.frame(cbind(Name=names(x), tbl))
  return(tbl)
}

allSlotsEqual <-
# For a list of objects and a slot, return logical indicating whether
# all slots equal value
function(x, slot, value=TRUE) {
  return(all(sapply(x, function(y) slot(y, slot))))
}

report <-
# report generates a report using knitr. It accepts a list of objects
# or a single object that inherit from SequenceSummary. In the latter
# case, we make a single list of it; this allows downstream functions
# to work on a consistent object class. Using a generic method is not
# really needed.
function(x, outdir=".", kmers=TRUE, hash=FALSE, landscape=TRUE) {
  if ((hash && !allSlotsEqual(x, "hashed")))
    stop("hash is TRUE but not all list elements have ")
  if ((hash && !allSlotsEqual(x, "kmers.hashed")))
    stop("if kmers or hash are TRUE")

  ll <- list()
  if (is(x, "SequenceSummary"))
    ll[[basename(x@filename)]] <- x
  else if (is.list(x))
    ll <- x
  else
    stop("report() only handles lists and objects that inherit from SequenceSummary")

  if (landscape) {
    message("using landscape layout")
    FIG.HEIGHT <- 7
    FIG.WIDTH <- 10
  }
  local({
    x <- ll
    include.kmers <- kmers
    include.hash <- hash
    knit(TEMPLATES$PDF)
  })
}

makeHashTable <- function(x, n=4) {
  tbl.list <- mapply(function(y, name) {
    if (!y@hashed)
      return()
    d <- y@hash[1:n]
    tbl <- data.frame(Name=name, Sequence=names(d),
                      Count=d, "Proportion of Total"=d/sum(y@seq.lengths),
                      check.names=FALSE)
    rownames(tbl) <- NULL
    tbl
  }, x, names(x), SIMPLIFY=FALSE)
  tbl <- do.call(rbind, tbl.list)
  tbl <- as.data.frame(tbl)
  return(tbl)
}
