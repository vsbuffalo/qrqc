## knitr-method.R -- new report generation using knitr
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

setMethod(report, "SequenceSummary",
          function(x, outdir=datedir(), type=c("tex", "html"),
                   kmers=TRUE, hash=FALSE, landscape=TRUE,
                   template=NULL, tar=FALSE, compression=c("none", "gzip", "bzip2"), quiet=FALSE) {
            x.listed <- list()
            x.listed[[basename(x@filename)]] <- x
            generateReport(x.listed, outdir, type, kmers, hash, landscape, template, tar, compression, quiet)
          })

setMethod(report, "list",
          function(x, outdir=datedir(), type=c("tex", "html"),
                   kmers=TRUE, hash=FALSE, landscape=TRUE,
                   template=NULL, tar=FALSE, compression=c("none", "gzip", "bzip2"), quiet=FALSE) {
            generateReport(x, outdir, type, kmers, hash, landscape, template, tar, compression, quiet)
          })


generateReport <-
# report generates a report using knitr. It accepts a list of objects
# or a single object that inherit from SequenceSummary. In the latter
# case, we make a single list of it; this allows downstream functions
# to work on a consistent object class. Using a generic method is not
# really needed.
function(x, outdir=datedir(), type=c("tex", "html"),
         kmers=TRUE, hash=FALSE, landscape=TRUE,
         template=NULL, tar=FALSE, compression=c("none", "gzip", "bzip2"),
         quiet=FALSE) {
  type <- match.arg(type)
  orientation <- c("portrait", "landscape")[as.integer(landscape) + 1]
  compression <- match.arg(compression)
    
  if ((hash && !allSlotsEqual(x, "hashed")))
    stop("hash is TRUE but not all list elements have ")
  if ((hash && !allSlotsEqual(x, "kmers.hashed")))
    stop("if kmers or hash are TRUE")

  if (compression != "none" && !tar)
    stop("compression argument only applicable if tar=TRUE")

  if (!landscape && type == "html")
    stop("landscap=TRUE is only available option for type='html'")

  outfile <- sprintf("qrqc-report.%s", type)

  if (!file.exists(outdir)) {
    message(sprintf("creating directory %s", outdir))
    dir.create(outdir)
  }
  
  local({
    include.kmers <- kmers
    include.hash <- hash

    oldwd <- getwd()
    setwd(outdir)
    if (is.null(template))
      template <- qrqc_options$templates[[type]]
    if (quiet)
      knitfun <- function(x, output) suppressMessages(capture.output(knit(x, output)))
    else
      knitfun <- knit
    knitfun(template, output=outfile)
    setwd(oldwd)
  })

  if (tar) {
    ext <- c(none=".tar", gzip=".tar.gz", bzip2=".tar.bz2")[compression]
    tarfile <- paste0(outdir, ext)
    message(sprintf("creating %s", tarfile))
    tar(tarfile, outdir, compression=compression)
    outfile <- tarfile
  }
  
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
