## ggplotting-methods.R -- new plotting functions that will coexist
## with existing qrqc plotting functions (however these others will
## warn that they are deprecated).

setMethod("getQual", signature(x="FASTQSummary"), 
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
})

setMethod("getGC", signature(x="SequenceSummary"), 
# Given a SequenceSummary object, extract a dataframe of GC content.
function(x) {
  gcd <- local({
    base.props <- getBaseProps(x)
    d <- subset(base.props, base.props$base %in% c('G', 'C'))
    gc <- aggregate(d$proportion, list(d$position), sum)
    colnames(gc) <- c('position', 'gc')
    return(gc)
  })
})

setMethod("getBase", signature(x="SequenceSummary"), 
# Given a SequenceSummary object, extract a dataframe of base content.
function(x, drop=TRUE) {
  base.freqs <- local({
    tmp <- melt(x@base.freqs, id='position')
    colnames(tmp) <- c('position', 'base', 'frequency')
    return(tmp)
  })
  if (drop) {
    # is zero for all bases?
    bt <- aggregate(base.freqs$frequency, list(base=base.freqs$base), sum)
    base.keep <- bt$base[bt$x > 0]
    return(subset(base.freqs, base %in% base.keep))
  }
  base.freqs
})

setMethod("getBaseProp", signature(x="SequenceSummary"),
# Given an object that inherits from SummarySequence, return the bases
# by proportion.
function(x, drop=TRUE) {
  base.props <- local({
    tmp <- prop.table(as.matrix(x@base.freqs[, -1]), margin=1)
    tmp <- melt(tmp, id='position')
    colnames(tmp) <- c('position', 'base', 'proportion')
    return(tmp)
  })
  if (drop) {
    # is zero for all bases?
    bt <- aggregate(base.props$proportion, list(base=base.props$base), sum)
    base.keep <- bt$base[bt$x > 0]
    return(subset(base.props, base %in% base.keep))
  }
  return(base.props)
})

setMethod("getSeqlen", signature(x="SequenceSummary"),
# Get sequence length dataframe from object that inherits from
# SequenceSummary.
function(x) {
  l <- x@seq.lengths
  x <- 1:length(l)
  data.frame(length=x, count=l)
})

setMethod("getMCQual", signature(x="FASTQSummary"),
# Generate random draws from a FASTQSummary object so that a lowess
# can be fit to quickly view the general trend of a quality summary.
function(x, n=100) {
  set.seed(0)
  vals <- as.numeric(colnames(x@qual.freqs)[-1])
  binsample <- function(x) {
    sample(vals, n, prob=x/sum(x), replace=TRUE)
  }

  # Use binsample() to sample from the binned quality frequency data, then
  # do some data munging to get it in a usable format.
  s <- apply(x@qual.freqs[, -1], 1, binsample)
  tmp <- t(s)
  nr <- nrow(tmp)
  dim(tmp) <- c(ncol(tmp)*nrow(tmp), 1)
  data.frame(position=1:nr, quality=tmp)
})

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
function(x, smooth=TRUE, extreme.color="grey", quantile.color="orange",
         mean.color="blue", median.color=NULL) {
  qd <- getQual(x)
  p <- ggplot(qd)
  p <- p + geom_qlinerange(extreme.color=extreme.color, quantile.color=quantile.color,
                           mean.color=mean.color, median.color=median.color)
  p <- p + scale_y_continuous("quality")
  if (smooth) {
    mcd <- getMCQual(x)
    p <- p + geom_smooth(aes(x=position, y=quality), data=mcd, se=FALSE)
  }
  p
})

setMethod("qualPlot", signature(x="list"),
# Plot a list of FASTQSummary objects as facets.
function(x, smooth=TRUE, extreme.color="grey", quantile.color="orange",
         mean.color="blue", median.color=NULL) {
  if (!length(names(x)))
    stop("A list pased into qualPlot must have named elements.")
  if (!all(sapply(x, class) == "FASTQSummary"))
    stop("All items in list must have class FASTQSummary.")
  qd <- list2df(x, getQual)
  p <- ggplot(qd)
  p <- p + geom_qlinerange(extreme.color=extreme.color, quantile.color=quantile.color,
                           mean.color=mean.color, median.color=median.color)
  p <- p + scale_y_continuous("quality")
  if (smooth) {
    mcd <- list2df(x, getMCQual)
    p <- p + geom_smooth(aes(x=position, y=quality), data=mcd, se=FALSE)
  }
  p <- p + facet_wrap( ~ sample)
  p
})


setMethod("gcPlot", signature(x="SequenceSummary"),
# Plot GC by read position.
function(x, color="red") {
  gcd <- getGC(x)
  p <- ggplot(gcd) + geom_line(aes(x=position, y=gc), color=color)
  p
})

setMethod("gcPlot", signature(x="list"),
# Plot GC by read position.
function(x, color="red") {
  gcd <- list2df(x, getGC)
  p <- ggplot(gcd) + geom_line(aes(x=position, y=gc), color=color) + facet_wrap( ~ sample)
  p
})

setMethod("basePlot", signature(x="SequenceSummary"),
# Plot Bases by read position.
function(x, geom=c("line", "bar", "dodge"), type=c("frequency", "proportion"), 
         bases=DNA_BASES_N, colorvalues=getBioColor("DNA_BASES_N")) {
  colorvalues <- colorvalues[bases]
  # get the type
  type <- match.arg(type)
  fun <- list(frequency=getBase, proportion=getBaseProp)[[type]]

  # get the geom
  geom <- match.arg(geom)
  geom.list <- list(line=geom_line(aes_string(x="position", y=type, color="base")),
                    bar=geom_bar(aes_string(x="position", y=type, fill="base"), stat="identity"),
                    dodge=geom_bar(aes_string(x="position", y=type, fill="base"), stat="identity", position="dodge"))
  g <- geom.list[[geom]]

  bd <- fun(x)

  p <- ggplot(subset(bd, base %in% bases)) + g
  p <- p + scale_colour_manual(values=colorvalues)
  p <- p + scale_fill_manual(values=colorvalues)
  p
})


setMethod("basePlot", signature(x="list"),
# Plot Bases by read position for a list.
function(x, geom=c("line", "bar", "dodge"), type=c("frequency", "proportion"),
         bases=DNA_BASES_N, colorvalues=getBioColor("DNA_BASES_N")) {
  colorvalues <- colorvalues[bases]
  # get the type
  type <- match.arg(type)
  fun <- list(frequency=getBase, proportion=getBaseProp)[[type]]
  bd <- list2df(x, fun)

  # get the geom
  geom <- match.arg(geom)
  geom.list <- list(line=geom_line(aes_string(x="position", y=type, color="base")),
                    bar=geom_bar(aes_string(x="position", y=type, fill="base"), stat="identity"),
                    dodge=geom_bar(aes_string(x="position", y=type, fill="base"), stat="identity", position="dodge"))
  g <- geom.list[[geom]]

  p <- ggplot(subset(bd, base %in% bases)) + g + facet_wrap( ~ sample)
  p <- p + scale_colour_manual(values=colorvalues)
  p <- p + scale_fill_manual(values=colorvalues)
  p
})

setMethod("seqlenPlot", signature(x="SequenceSummary"),
# Plot sequence lengths.
function(x) {
  ld <- getSeqlen(x)
  p <- ggplot(ld, aes(x=length, y=count)) + geom_bar(stat="identity")
  p
})


setMethod("seqlenPlot", signature(x="list"),
# Plot sequence lengths.
function(x) {
  ld <- list2df(x, getSeqlen)
  p <- ggplot(ld, aes(x=length, y=count)) + geom_bar(stat="identity")
  p <- p + facet_wrap( ~ sample)
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
