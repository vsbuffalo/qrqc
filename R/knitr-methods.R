## knitr-method.R -- new report generation using knitr
TEMPLATES <- list(tex=system.file('extdata', 'knitr-report-template.Rnw', package='qrqc'),
                  html=system.file('extdata', 'knitr-report-template.html', package='qrqc'))
FIG.HEIGHT <- 6
FIG.WIDTH <- 8

seqLengthRange <- function(x) {
  s <- which(x@seq.lengths > 0)
  return(c(min(s), max(s)))
}

reportName <- function(x) {
  if (is.list(x)) {
    tmp <- sapply(x, function(y) strsplit(basename(y@filename), '.', fixed=TRUE)[[1]][1])
    out <- paste(tmp, sep="-")
  } else {
    out <- basename(y@filename)
  }
  return(out)
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

datedir <-
# make an output directory based on the current time. Checks if
# directory exists and errors out if it would overwrite it.
function(prefix="qrqc-report-") {
  dirname <- paste0(prefix, format(Sys.time(), "%Y-%m-%dT%H%M%S"))
  if (file.exists(dirname))
    stop(sprintf("directory '%s' already exists!", dirname))
  dir.create(dirname)
  return(dirname)
}

setMethod(report, "SequenceSummary",
          function(x, outdir=datedir(), type=c("tex", "html"),
                   kmers=TRUE, hash=FALSE, landscape=TRUE,
                   template=NULL) {
            x.listed <- list()
            x.listed[[basename(x@filename)]] <- x
            generateReport(x.listed, outdir, type, kmers, hash, landscape, template)
          })

setMethod(report, "list",
          function(x, outdir=datedir(), type=c("tex", "html"),
                   kmers=TRUE, hash=FALSE, landscape=TRUE,
                   template=NULL) {
            generateReport(x, outdir, type, kmers, hash, landscape, template)
          })


generateReport <-
# report generates a report using knitr. It accepts a list of objects
# or a single object that inherit from SequenceSummary. In the latter
# case, we make a single list of it; this allows downstream functions
# to work on a consistent object class. Using a generic method is not
# really needed.
function(x, outdir=datedir(), type=c("tex", "html"),
         kmers=TRUE, hash=FALSE, landscape=TRUE,
         template=NULL) {
  type <- match.arg(type)

  if ((hash && !allSlotsEqual(x, "hashed")))
    stop("hash is TRUE but not all list elements have ")
  if ((hash && !allSlotsEqual(x, "kmers.hashed")))
    stop("if kmers or hash are TRUE")

  if (landscape) {
    FIG.HEIGHT <- 7
    FIG.WIDTH <- 10
  }

  outfile <- file.path(outdir, sprintf("qrqc-report.%s", type))

  local({
    include.kmers <- kmers
    include.hash <- hash

    if (is.null(template))
      template <- TEMPLATES[[type]]
    knit(template, output=outfile)
  })
  return(outfile)
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
