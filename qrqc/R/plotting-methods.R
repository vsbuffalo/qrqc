## plotting-methods.R - functions to create plots from FASTQSummary and FASTASummary objects

setMethod("plotBases", "SequenceSummary",
# Plot the frequency (absolute counts) and proportion of bases across
# a read.
function(obj, type="freq", bases=NULL, legend=TRUE) {
  if (!(type %in% c("freq", "prop")))
    stop("'type' must be either 'freq' or 'prop'.")

  plot.new()

  if (type == "freq") {
    base.freqs <- local({
      tmp <- melt(obj@base.freqs, id='position')
      colnames(tmp) <- c('position', 'base', 'frequency')
      return(tmp)
    })

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

    if (legend) {
      # add legend
      legend.bases <- levels(base.freqs$base[which(base.freqs$frequency > 0), drop=TRUE])
      legend.colors <- NUCLEOTIDES.COLORS[as.character(legend.bases)]
      legend("topright",
             legend.bases, col=legend.colors, lwd=2, box.lwd=0.6)
    }
    
  } else if (type == "prop") {
    base.props <- getBaseProps(obj)
    
    w <- 1.2
    plot.window(ylim=c(min(base.props$proportion)*w, max(base.props$proportion)*w),
                xlim=c(1, max(base.props$position)))
    
    if (is.null(bases))
      levs <- levels(base.props$base)
    else
      levs = bases

    for (i in 1:length(levs)) {
      base <- levs[i]
      lines(base.props$position[base.props$base == base],
            base.props$proportion[base.props$base == base],
            col=NUCLEOTIDES.COLORS[as.character(base)])
      }
    
    axis(1, at=min(base.props$position):max(base.props$position))
    axis(2, at=seq(0, 1, by=0.2))
    abline(h=0.25, col='grey')  
    title(main="base proportion by position in read", xlab="position", ylab="proportion")

    if (legend) {
      # add legend
      legend.bases <- levels(base.props$base[which(base.props$proportion > 0), drop=TRUE])
      legend.colors <- NUCLEOTIDES.COLORS[as.character(legend.bases)]
      legend("topright",
             legend.bases, col=legend.colors, lwd=2, box.lwd=0.6)
    }
  }
  
})


setMethod("plotQuals", "FASTQSummary",
# Make a series of (custom) boxplots, with a point for median and a
# horizontal line for the mean. If `lowess` is TRUE, add lowess curve,
# which is fit through MC samples though the binned quals with
# `qualMCLowess`.
function(obj, ylim='relative', lowess=TRUE, histogram=TRUE, legend=FALSE) {
  d <- local({
    tmp <- apply(obj@qual.freqs[, -1], 1, binned2boxplot)
    tmp <- t(tmp)
    tmp <- cbind(1:nrow(tmp), tmp)
    vals <- as.numeric(names(obj@qual.freqs[, -1]))
    tmp <- cbind(tmp, apply(obj@qual.freqs[, -1], 1, function(x) weighted.mean(vals, x)))
    tmp <- as.data.frame(tmp)
    colnames(tmp) <- c('position', 'ymin', 'alt.lower', 'lower',
                       'middle', 'upper', 'alt.upper', 'ymax', 'mean')
    return(tmp)
  })

  if (ylim == 'fixed') {
    q <- QUALITY.CONSTANTS[[obj@quality]]
    qmin <- q$min
    qmax <- q$max
  } else {
    qmin <- min(d$ymin)
    qmax <- max(d$ymax)
  }

  op <- par(no.readonly=TRUE)
  m <- sprintf("quality distribution by read base\n(quality type: %s)", obj@quality)
  
  if (histogram) {
    nf <- layout(c(1, 2), heights=c(1, 4))
    #layout.show(nf)
      
    ## First plot: histogram of sequence lengths
    par(mar=c(0, 4, 3, 1))
    s <- obj@seq.lengths
    plot.new()
    plot.window(xlim=c(0, nrow(d)), ylim=c(0, 1))
    axis(1, at=1:nrow(d))
    axis(2, at=seq(0, 1, by=0.1))
    lines(prop.table(s[1:max(which(s != 0))]), type='h', col='blue', lwd=2)
    m <- sprintf("quality distribution by read base and sequence length histogram\n(quality type: %s)", obj@quality)
    title(main=m, ylab="density")
    par(mar=c(4.5, 4, 1, 1)) # setup for next plot
  }

  ## Second plot: quantile plots  
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
  if (!histogram)
    title(main=m)
  par(op)

  if (legend) {
    legend("bottomleft",
           c("10% - 90%", "25% - 75%", "median", "mean", "lowess"),
           col=c("grey", "orange", "blue", "dark green", "purple"), lwd=2, box.lwd=0.6)
  }
})

setMethod("plotGC", "SequenceSummary",
# Plot the proportion of bases that are G or C (the GC content).
function(obj) {
  gc <- local({
    base.props <- getBaseProps(obj)
    d <- subset(base.props, base.props$base %in% c('G', 'C'))
    gc <- aggregate(d$proportion, list(d$position), sum)
    colnames(gc) <- c('position', 'gc')
    return(gc)
  })
  
  plot.new()

  w <- 1.2
  plot.window(ylim=c(min(gc$gc)*w, max(gc$gc)*w),
              xlim=c(1, max(gc$position)))
  
  
  lines(gc$position, gc$gc, col='pink')

  axis(1, at=min(gc$position):max(gc$position))
  axis(2, at=seq(0, 1, by=0.2))
  abline(h=0.5, col='grey')  
  title(main="gc content by position in read", xlab="position", ylab="proportion")
  
})

setMethod("plotSeqLengths", "SequenceSummary",
# Plot a histogram of sequence lengths. This is done manually (not
# using hist) so that the case when there's a single length, the
# histogram doesn't look ridiculous.
function(obj) {
  plot.new()

  l <- obj@seq.lengths
  x <- unique(1:length(l))
  d <- data.frame(cbind(x, l))
  d$p <- d$l/sum(d$l)
  browser()
  plot.window(ylim=c(0, 1), xlim=c(min(x), max(x)))
  apply(d, 1, function(x) lines(x[1], x[3], type='h', col='blue', lwd=2))
  axis(1, at=(min(x)-1):(max(x)+1))

  axis(2, at=seq(0, max(x), by=0.2))
  title(ylab='density', xlab='sequence length',
        main='sequence length distribution')
})
