## utils.R - Utility functions for creating parts of graphics, loading data, etc.

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
function(m) {
  cs <- colSums(m)
  return(m[, 1:max(which(cs != 0))])
}

.trimArray <-
function(x)
  x[1:max(which(x != 0))]

sortSequenceHash <- function(seq.hash) {
  return(sort(unlist(seq.hash), decreasing=TRUE))
}

lengths2weights <-
# Given a histogram of seq.lengths, create a weighting scheme for
# estimating mean quality from sequences of different lengths.
# This converts a histogram to a base-by-base coverage vector.
function(length.dist) {
  last <- length(.trimArray(length.dist))
  l <- length.dist[1:last]
  s <- numeric(last)
  for (i in last:1) {
    s[1:i] <- s[1:i] + l[i]
  }
  return(s)
}

meanFromBins <- 
# Given binned quality data, return the mean quality.
function(qual.dist, length.dist=NULL) {
  means <- numeric(nrow(qual.dist))
  for (i in 1:nrow(qual.dist)) {
    vals <- as.integer(colnames(qual.dist[-1]))
    means[i] <- weighted.mean(vals, w=qual.dist[i, -1])
  }
  if (is.null(length.dist))
    return(mean(means))
  return(weighted.mean(means, w=lengths2weights(length.dist)))
}

binned2quantilefunc <-
# Assumes names are x values
function(x) {
  y <- as.numeric(x[x!=0])
  x <- as.integer(names(x)[x!=0])
  if (length(x) < 10) {
    # too few values to interpolate ECDF
    return(NULL)
  }
  t <- approxfun(cumsum(y)/sum(y), x, yleft=min(x), yright=max(x))
  return(t)
}

binned2boxplot <-
# Extract useful quantiles for a boxplot from an empirical qunatile
# function.
function(x) {
  f <- binned2quantilefunc(x)
  if (is.null(f)) {
    # Too few values to interpolate; build actual vector (since it's
    # small enough to do so) and find certain percentiles.
    y <- as.numeric(x[x!=0])
    x <- as.integer(names(x)[x!=0])
    d <- quantile(rep(x, times=y), probs=c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1))
    out <- c(ymin=d[1],
             alt.lower=d[2],
             lower=d[3],
             middle=d[4],
             upper=d[5],
             alt.upper=d[6],
             ymax=d[7])
    return(out)
  }

  # binned2quantilefunc had enough values to interpolate ECDF.
  out <- c(ymin=f(0),
           alt.lower=f(0.1),
           lower=f(0.25),
           middle=f(0.5),
           upper=f(0.75),
           alt.upper=f(0.9),
           ymax=f(1))
  return(out)
}

qualMCLowess <-
# Monte Carlo Lowess, or lowess on binned quality data
function(obj, n=100, f=1/6) {
  set.seed(0)
  vals <- as.numeric(names(obj@qual.freqs[1, -1]))
  binsample <- function(x) {
    sample(vals, n, prob=x/sum(x), replace=TRUE)
  }

  # Use binsample() to sample from the binned quality frequency data, then
  # do some data munging to get it in a usable format for lowess.
  d <- local({
    s <- apply(obj@qual.freqs[, -1], 1, binsample)
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

getBaseProps <-
# Given an object that inherits from SummarySequence, return the bases
# by proportion.
function(obj) {
  base.props <- local({
    tmp <- prop.table(as.matrix(obj@base.freqs[, -1]), margin=1)
    tmp <- melt(tmp, id='position')
    colnames(tmp) <- c('position', 'base', 'proportion')
    return(tmp)
  })
  return(base.props)
}
