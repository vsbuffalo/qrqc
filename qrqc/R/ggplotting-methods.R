## ggplotting-methods.R -- new plotting functions that will coexist
## with existing qrqc plotting functions (however these others will
## warn that they are deprecated).

## Methods for plotting FASTQSummary, FASTASummary, and
## SequenceSummary objects


getQual <- 
# Given a FASTQSummary object, extract a dataframe of quality
# distribution statistics using binned2boxplot.
function(x) {
  tmp <- apply(x@qual.freqs[, -1], 1, binned2boxplot)
  tmp <- t(tmp)
  tmp <- cbind(1:nrow(tmp), tmp)
  vals <- as.numeric(names(x@qual.freqs[, -1]))
  tmp <- cbind(tmp, apply(x@qual.freqs[, -1], 1, function(x) weighted.mean(vals, x)))
  tmp <- as.data.frame(tmp)
  colnames(tmp) <- c('position', 'ymin', 'alt.lower', 'lower',
                     'middle', 'upper', 'alt.upper', 'ymax', 'mean')
  return(tmp)
}

getGC <-
# Give a SequenceSummary object, extract a dataframe of GC content
function(x) {
  gcd <- local({
    base.props <- getBaseProps(x)
    d <- subset(base.props, base.props$base %in% c('G', 'C'))
    gc <- aggregate(d$proportion, list(d$position), sum)
    colnames(gc) <- c('position', 'gc')
    return(gc)
  })
}

list2df <-
# Use mapply
function(x, fun) {
  d <- mapply(function(x, n) {
    d.tmp <- fun(x)
    d.tmp$sample <- n
    d.tmp
  }, x, names(x), SIMPLIFY=FALSE)
  as.data.frame(do.call(rbind, d))
}
  
setMethod("qualPlot", signature(x="FASTQSummary"),
# Plot a single FASTQSummary object.
function(x, extreme.color="grey", quantile.color="orange",
         mean.color="blue", median.color=NULL, ...) {
  qd <- getQual(x)
  p <- ggplot(qd)
  p <- p + geom_qlinerange(extreme.color=extreme.color, quantile.color=quantile.color,
                           mean.color=mean.color, median.color=median.color)
  p <- p + scale_y_continuous("quality")
  p
})

setMethod("qualPlot", signature(x="list"),
# Plot a list of FASTQSummary objects as facets.
function(x, extreme.color="grey", quantile.color="orange",
         mean.color="blue", median.color=NULL, ...) {
  if (!length(names(x)))
    stop("A list pased into qualPlot must have named elements.")
  qd <- list2df(x, getQual)
  p <- ggplot(qd)
  p <- p + geom_qlinerange(extreme.color=extreme.color, quantile.color=quantile.color,
                           mean.color=mean.color, median.color=median.color)
  p <- p + facet_grid(. ~ sample) + scale_y_continuous("quality")
  p
})


setMethod("gcPlot", signature=(x="SequenceSummary"),
# Plot GC by base
function(x, color="red") {
  gcd <- getGC(x)
  p <- ggplot(gcd) + geom_line(aes(x=position, y=gc), color=color)
  p
})

setMethod("gcPlot", signature=(x="list"),
# Plot GC by base
function(x, color="red") {
  gcd <- list2df(x, getGC)
  p <- ggplot(gcd) + geom_line(aes(x=position, y=gc), color=color) + facet_grid(. ~ sample)
  p
})

geom_qlinerange <- 
# Add a series of geoms to plot quality statistics.
function(extreme.color="grey", quantile.color="orange", mean.color="blue", median.color=NULL) {
  args <- as.list(match.call(call = sys.call(sys.parent()))[-1])
  l <- list(extreme.color=geom_linerange(aes(x=position, ymin=alt.lower, ymax=alt.upper), color=extreme.color),
            quantile.color=geom_linerange(aes(x=position, ymin=lower, ymax=upper, y=mean), color=quantile.color, size=1.2),
            mean.color=geom_point(aes(x=position, y=mean), color=mean.color),
            median.color=geom_point(aes(x=position, y=middle), color=median.color))

  keep <- sapply(names(l), function(n) !is.null(get(n)))
  l[keep]
}





