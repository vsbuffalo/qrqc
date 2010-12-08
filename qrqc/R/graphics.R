## graphics.R - functions to create plots from SeqStats objects

plotBaseFreqs <-
# Plot the frequency (absolute counts) of bases across a read.
function(obj, bases=NULL) { 
  base.freqs <- local({
    tmp <- melt(obj@base.freqs, id='position')
    colnames(tmp) <- c('position', 'base', 'frequency')
    return(tmp)
  })

  plot.new()
  plot.window(ylim=c(min(base.freqs$frequency)*1.2, max(base.freqs$frequency)*1.2),
              xlim=c(1, max(base.freqs$position)))

  if (is.null(bases))
    levs <- levels(base.freqs$base)
  else
    levs = bases

  for (base in levs) {
    lines(base.freqs$position[base.freqs$base == base],
          base.freqs$frequency[base.freqs$base == base],
          col=NUCLEOTIDES.COLORS[as.character(base)])
  }

  axis(1, at=min(base.freqs$position):max(base.freqs$position))
  axis(2)
  title(main="base frequency by position in read", xlab="position", ylab="frequency")
}

plotBaseProps <-
# Plot proportions of all nucleotides across a read.
function(obj) {
  base.props <- obj@base.props
  plot.new()

  w <- 1.2
  plot.window(ylim=c(min(base.props$proportion)*w, max(base.props$proportion)*w),
              xlim=c(1, max(base.props$position)))
  
  
  for (base in levels(base.props$base)) {
    lines(base.props$position[base.props$base == base],
          base.props$proportion[base.props$base == base],
          col=NUCLEOTIDES.COLORS[as.character(base)])
  }

  axis(1, at=min(base.props$position):max(base.props$position))
  axis(2, at=seq(0, 1, by=0.2))
  abline(h=0.25, col='grey')  
  title(main="base proportion by position in read", xlab="position", ylab="proportion")
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
# Extract useful quantiles for a boxplot from an empirical qunatile
# function.
function(x) {
  f <- binned2quantilefunc(x)
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


plotQuals <-
# Make a series of (custom) boxplots, with a point for median and a
# horizontal line for the mean. If `lowess` is TRUE, add lowess curve,
# which is fit through MC samples though the binned quals with
# `qualMCLowess`.
function(obj, ylim='relative', lowess=TRUE) {
  if (obj@type == 'fasta')
    stop("plotQuals only works on SeqStat objects with quality (i.e. from a FASTA file)")
  d <- local({
    tmp <- apply(obj@qual.freqs[, -1], 1, binned2boxplot)
    tmp <- t(tmp)
    tmp <- cbind(1:nrow(tmp), tmp)
    vals <- as.numeric(names(obj@qual.freqs[, -1]))
    tmp <- cbind(tmp, apply(obj@qual.freqs[, -1], 1, function(x) weighted.mean(vals, x)))
    colnames(tmp)[1] <- 'position'
    colnames(tmp)[9] <- 'means'
    return(as.data.frame(tmp))
  })

  if (ylim == 'fixed') {
    q <- QUALITY.CONSTANTS[[obj@quality]]
    qmin <- q$min
    qmax <- q$max
  } else {
    qmin <- min(d$ymin)
    qmax <- min(d$ymax)
  }

  op <- par(no.readonly=TRUE)
  nf <- layout(c(1, 2), heights=c(1, 4))
  #layout.show(nf)

  ## First plot: histogram of sequence lengths
  par(mar=c(0, 4, 3, 1))
  s <- obj@seq.lengths
  plot.new()
  plot.window(xlim=c(0, nrow(d)), ylim=c(0, 1))
  axis(1, at=1:nrow(d))
  axis(2, at=seq(0, 1, by=0.1))
  lines(prop.table(s[2:max(which(s != 0))]), type='h', col='blue', lwd=2) # offset by one, since C uses 0-indexin
  m <- sprintf("quality distribution by read base and sequence length histogram (quality type: %s)", obj@quality)
  title(main=m, ylab='density')

  ## Second plot: quantile plots
  par(mar=c(4.5, 4, 1, 1))
  plot.new()
  plot.window(xlim=c(0, nrow(d)), ylim=c(qmin, qmax))

  axis(1, at=1:nrow(d), las=1)
  axis(2, at=qmin:qmax, las=1)

  apply(d, 1, function(x)
        lines(x=c(x[1], x[1]), y=c(x[3], x[7]), col='grey'))
  
  apply(d, 1, function(x)
        lines(x=c(x[1], x[1]), y=c(x[4], x[6]), lwd=2.5, col='orange'))
  points(d$position, d$middle, pch=20, col='blue')

  lw <- 0.2
  apply(d, 1, function(x)
        lines(x=c(x[1]-lw, x[1]+lw), y=c(x[9], x[9]), lwd=2.5, col='dark green'))

  if (lowess)
    lines(qualMCLowess(obj), col='purple')

  title(xlab="position", ylab="quality")
  par(op)
}

plotGC <-
# Plot the proportion of bases that are G or C (the GC content).
function(obj) {
  gc <- local({
    d <- subset(obj@base.props, obj@base.props$base %in% c('G', 'C'))
    gc <- aggregate(d$proportion, list(d$position), sum)
    colnames(gc) <- c('position', 'gc')
    return(gc)
  })
  
  plot.new()

  w <- 1.2
  plot.window(ylim=c(min(gc$gc)*w, max(gc$gc)*w),
              xlim=c(1, max(gc$position)))
  
  
  lines(gc$position, gc$gc, col='red')

  axis(1, at=min(gc$position):max(gc$position))
  axis(2, at=seq(0, 1, by=0.2))
  abline(h=0.5, col='grey')  
  title(main="gc content by position in read", xlab="position", ylab="proportion")
  
}

plotSeqLengths <-
# Plot a histogram of sequence lengths. This is done manually (not
# using hist) so that the case when there's a single length, the
# histogram doesn't look ridiculous.
function(obj) {
  plot.new()

  l <- obj@seq.lengths
  x <- unique(1:length(l))
  d <- data.frame(cbind(x, l))
  d$p <- d$l/sum(d$l)

  plot.window(ylim=c(0, 1), xlim=c(min(x), max(x)))
  apply(d, 1, function(x) lines(x[1], x[2], type='h', col='blue', lwd=2))
  axis(1, at=(min(x)-1):(max(x)+1))

  axis(2, at=seq(0, max(x), by=0.2))
  title(ylab='density', xlab='sequence length',
        main='sequence length distribution')
}
