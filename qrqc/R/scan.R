## scan.R - Scan sequenece files

QUALITY.CONSTANTS <- list(phred=list(offset=33, min=0, max=93),
                          solexa=list(offset=64, min=-5, max=62),
                          illumina=list(offset=64, min=0, max=62))

NUCLEOTIDES <- c('A', 'T', 'C', 'G', 'N')

.setQualityNames <-
# Given a quality type (as integer), name the matrix output from
# the C function `summarize_fastq_file` accordingly.
function(matrix, quality) {
  constants <- QUALITY.CONSTANTS[[quality]]
  rownames(matrix) <- constants$min:constants$max
  return(matrix)
}

.trimRightCols <-
# Remove blank cols from right of matrix, which occur because the
# matrix allocates excess space.
function(matrix) {
  cs <- colSums(matrix)
  return(matrix[, 1:max(which(cs != 0))])
}

sortSequenceHash <- function(seq.hash) {
  return(sort(unlist(seq.hash), decreasing=TRUE))
}

summarizeFastq <-
# Use the C function summarize_fastq_file to create matrices of base
# and quality counts, per position along all reads.
# TODO: Currently only the illumina quality setting has been tested.
function(filename, max.length=100, quality='illumina', hash=TRUE, verbose=FALSE) {
  if (!file.exists(filename))
    stop(sprintf("file '%s' does not exist", filename))

  out <- .Call('summarize_fastq_file', filename,
               as.integer(max.length),
               as.integer(which(names(QUALITY.CONSTANTS) == quality) - 1),
               as.logical(hash),
               as.logical(verbose))

  names(out) <- c('base.freqs', 'qual.freqs')

  if (hash) {
    names(out)[3] <- 'hash'
    out$hash <- sortSequenceHash(out$hash)
  }
  
  ## Data cleaning
  out$base.freqs <- local({
    tmp <- as.data.frame(t(.trimRightCols(out$base.freqs)))
    tmp <- cbind(1:nrow(tmp), tmp)
    colnames(tmp) <- c('position', NUCLEOTIDES)
    return(tmp)})

  out$qual.freqs <- local({
    tmp <- .trimRightCols(out$qual.freqs)
    tmp <- as.data.frame(t(.setQualityNames(tmp, quality)))
    tmp <- cbind(1:nrow(tmp), tmp)
    colnames(tmp)[1] <- 'position'
    return(tmp)})

  out$base.props <- local({
    tmp <- prop.table(as.matrix(out$base.freqs[, -1]), margin=1)
    tmp <- melt(tmp, id='position')
    colnames(tmp) <- c('position', 'base', 'proportion')
    return(tmp)
  })

  out$quality=quality
  return(out)
}

plotBaseFreqs <-
#
function(obj){
  base.freqs <- local({
    tmp <- melt(obj$base.freqs, id='position')
    colnames(tmp) <- c('position', 'base', 'frequency')
    return(tmp)
  })

  p <- ggplot(base.freqs, aes(x=position, y=frequency, color=base))
  p <- p + geom_line()
  print(p)
  invisible(p)
}

plotBaseProps <-
#
function(obj){
  p <- ggplot(obj$base.props, aes(x=position, y=proportion, color=base))
  p <- p + geom_line()
  print(p)
  invisible(p)  
}

binned2quantilefunc <-
# Assumes names are x values
function(x) {
  y <- as.numeric(x[x!=0])
  x <- as.integer(names(x)[x!=0])
  t <- approxfun(cumsum(y)/sum(y), x, yleft=min(x), yright=max(x))
  return(t)
}

binned2boxplot <-
#
function(x) {
  f <- binned2quantilefunc(x)
  out <- c(ymin=f(0),
           lower=f(0.25),
           middle=f(0.5),
           upper=f(0.75),
           ymax=f(1))
  return(out)
}

qualMCLowess <-
# Monte Carlo Lowess, or lowess on binned quality data
function(obj, n=100, f=1/6) {
  vals <- as.numeric(names(obj$qual.freqs[1, -1]))
  binsample <- function(x) {
    sample(vals, n, prob=x/sum(x), replace=TRUE)
  }

  # Use binsample() to sample from the binned quality frequency data, then
  # do some data munging to get it in a usable format for lowess.
  d <- local({
    s <- apply(obj$qual.freqs[, -1], 1, binsample)
    tmp <- t(s)
    nr <- nrow(tmp)
    dim(tmp) <- c(ncol(tmp)*nrow(tmp), 1)
    tmp <- as.data.frame(cbind(1:nr, tmp))    
    colnames(tmp) <- c('position', 'quality')
    return(tmp)
  })
  l <- lowess(d$position, d$quality, f=f)
  return(l)
}


plotQuals <-
#
function(obj, ylim='relative') {
  d <- local({
    tmp <- apply(obj$qual.freqs[, -1], 1, binned2boxplot)
    tmp <- t(tmp)
    tmp <- cbind(1:nrow(tmp), tmp)
    vals <- as.numeric(names(obj$qual.freqs[, -1]))
    tmp <- cbind(tmp, apply(obj$qual.freqs[, -1], 1, function(x) weighted.mean(vals, x)))
    colnames(tmp)[1] <- 'position'
    colnames(tmp)[7] <- 'means'
    return(as.data.frame(tmp))
  })

  if (ylim == 'fixed') {
    q <- QUALITY.CONSTANTS[[obj$quality]]
    qmin <- q$min
    qmax <- q$max
  } else {
    qmin <- min(d$ymin)
    qmax <- min(d$ymax)
  }
  
  plot.new()
  plot.window(xlim=c(0, nrow(d)), ylim=c(qmin, qmax))

  axis(1, at=1:nrow(d), las=1)
  axis(2, at=qmin:qmax, las=1)

  apply(d, 1, function(x) lines(x=c(x[1], x[1]), y=c(x[2], x[6]), col='grey'))
  apply(d, 1, function(x) lines(x=c(x[1], x[1]), y=c(x[3], x[5]), lwd=2.5, col='orange'))
  points(d$position, d$middle, pch=20, col='blue')

  lw <- 0.2
  apply(d, 1, function(x) lines(x=c(x[1]-lw, x[1]+lw), y=c(x[7], x[7]), lwd=2.5, col='green'))
}
plotQuals(s)
