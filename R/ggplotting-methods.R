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

setMethod("getKmer", signature(x="SequenceSummary"), 
# Given a SequenceSummary object, extract the k-mer dataframe
function(x) {
  x@kmer
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
  gcd
})

setMethod("getBase", signature(x="SequenceSummary"), 
# Given a SequenceSummary object, extract a dataframe of base content.
function(x, drop=TRUE) {
  base.freqs <- local({
    tmp <- reshape::melt(x@base.freqs, id='position')
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
    tmp <- reshape::melt(tmp, id='position')
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
  if (!length(names(x)))
    stop("list 'x' must have named elements.")
  d <- mapply(function(x, n) {
    d.tmp <- fun(x)
    d.tmp$sample <- n
    d.tmp
  }, x, names(x), SIMPLIFY=FALSE)
  as.data.frame(do.call(rbind, d))
}
  
setMethod("qualPlot", signature(x="FASTQSummary"),
# Plot a single FASTQSummary object.
function(x, smooth=TRUE, extreme.color="grey", quartile.color="orange",
         mean.color="blue", median.color=NULL) {
  qd <- getQual(x)
  p <- ggplot(qd)
  p <- p + geom_qlinerange(extreme.color=extreme.color, quartile.color=quartile.color,
                           mean.color=mean.color, median.color=median.color)
  p <- p + scale_y_continuous("quality")
  if (smooth) {
    mcd <- getMCQual(x)
    p <- p + geom_smooth(aes(x=position, y=quality), data=mcd, se=FALSE, color="purple")
  }
  p
})

setMethod("qualPlot", signature(x="list"),
# Plot a list of FASTQSummary objects as facets.
function(x, smooth=TRUE, extreme.color="grey", quartile.color="orange",
         mean.color="blue", median.color=NULL) {
  if (!length(names(x)))
    stop("A list pased into qualPlot must have named elements.")
  if (!all(sapply(x, class) == "FASTQSummary"))
    stop("All items in list must have class FASTQSummary.")
  qd <- list2df(x, getQual)
  p <- ggplot(qd)
  p <- p + geom_qlinerange(extreme.color=extreme.color, quartile.color=quartile.color,
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

  bd <- fun(x, drop=FALSE) # drop is FALSE if a user asks just to see N, and there are none

  p <- ggplot(subset(bd, base %in% bases)) + g

  if (geom %in% c("bar", "dodge"))
    p <- p + scale_fill_manual(values=colorvalues)
  if (geom == "line")
    p <- p + scale_color_manual(values=colorvalues)
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
  if (geom %in% c("bar", "dodge"))
    p <- p + scale_fill_manual(values=colorvalues)
  if (geom == "line")
    p <- p + scale_color_manual(values=colorvalues)
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
function(extreme.color="grey", quartile.color="orange", mean.color="blue", median.color=NULL) {
  l <- list(extreme.color=geom_linerange(aes_string(x="position", ymin="alt.lower", ymax="alt.upper"), color=extreme.color),
            quartile.color=geom_linerange(aes_string(x="position", ymin="lower", ymax="upper", y="mean"), color=quartile.color, size=1.2),
            mean.color=geom_point(aes_string(x="position", y="mean"), color=mean.color),
            median.color=geom_point(aes_string(x="position", y="middle"), color=median.color))

  keep <- sapply(names(l), function(n) !is.null(get(n)))
  l[keep]
}

scale_color_dna <-
# Make a color scheme using biovizBase's good choices.
function() {
  scale_color_manual(values=getBioColor("DNA_BASES_N"))
}

scale_color_iupac <-
# Make a color scheme using biovizBase's good choices.
function() {
  scale_color_manual(values=getBioColor("IUPAC_CODE_MAP"))
}

calcKL <- function(x) {
  d <- getKmer(x)

  kmerDist <- function(x) {
    # Given x, find the distribution of k-mers averaged across
    # position.
    tmp <- aggregate(x$count, list(kmer=x$kmer), sum)
    data.frame(kmer=tmp$kmer, prop=tmp$x/sum(tmp$x))    
  }
  
  kmer.dist <- kmerDist(d)
  kmer.by.pos <- split(d, list(position=d$position))

  ## We want to calculate the K-L divergence, but our sample spaces
  ## are not the same, as our k-mer for all positions is populated
  ## with more k-mers than most positions (as there aren't all k-mers
  ## in every position). Thus we calculate the probabilities of k-mers
  ## for all positions for the sub-sample space of the positional
  ## k-mers.
  y <- lapply(kmer.by.pos, kmerDist)
  
  kl <- mapply(function(x, p) {
    x$q <- kmer.dist$prop[match(x$kmer, kmer.dist$kmer)]
    klout <- data.frame(kmer=x$kmer, position=p,
                        kl=x$prop*log2(x$prop/x$q), p=x$prop, q=x$q)
    # by definition of K-L divergence, we only keep if for i, p(i) > 0
    # AND q(i) > 0.
    subset(klout, p > 0 & q > 0)
  }, y, as.numeric(names(y)), SIMPLIFY=FALSE)

  kld <- do.call(rbind, kl)
  kld
}


setMethod("kmerKLPlot", signature(x="SequenceSummary"),
# K-L divergence of empirical distribution across all reads to
# positional distribution.
function(x, n.kmers=20) {
  if (!nrow(getKmer(x)))
    stop("Data frame of k-mer counts by position is empty. Rerun readSeqFile with kmer=TRUE.")
  kld <- calcKL(x)  
  kld.subset <- subset(kld, kmer %in% kld$kmer[order(kld$kl, decreasing=TRUE)[seq_len(n.kmers)]])

  p <- ggplot(kld.subset)
  p <- p + geom_bar(aes(x=position, y=kl, fill=kmer), stat="identity")
  p + scale_y_continuous("K-L divergence of selected components")
})


setMethod("kmerKLPlot", signature(x="list"),
# K-L divergence of empirical distribution across all reads to
# positional distribution, facted by object.
function(x, n.kmers=20) {
  if (!all(sapply(x, function(x) nrow(getKmer(x)))))
    stop("Data frame of k-mer counts by position is empty. Rerun readSeqFile with kmer=TRUE.")

  kld.subset  <- list2df(x, function(x) {
    kld <- calcKL(x)  
    subset(kld, kmer %in% kld$kmer[order(kld$kl, decreasing=TRUE)[seq_len(n.kmers)]])
  })
  
  p <- ggplot(kld.subset)
  p <- p + geom_bar(aes(x=position, y=kl, fill=kmer), stat="identity")
  p + facet_wrap(~ sample) + scale_y_continuous("K-L divergence of selected components")
})


kmerEntropy <-
# Get a dataframe of k-mer entropy by position.
function(x) {
  d <- getKmer(x)

  y <- split(d, list(position=d$position))
  kmer.entropies <- sapply(y, function(x) {
    p <- x$count/sum(x$count)
    -sum(p*log2(p))
  })
  data.frame(position=as.numeric(names(y)), entropy=kmer.entropies)
}

setMethod("kmerEntropyPlot", signature(x="SequenceSummary"),
# Plot the k-mer entropy by position.
function(x) {
  ggplot(data=kmerEntropy(x)) + geom_line(aes(x=position, y=entropy), color="blue")
})

setMethod("kmerEntropyPlot", signature(x="list"),
# Plot the k-mer entropy by position.
function(x) {
  ke <- list2df(x, kmerEntropy)
  p <- ggplot(data=ke) + geom_line(aes(x=position, y=entropy), color="blue")
  p <- p + facet_wrap(~ sample)
  p
})

