## main.R - functions to handle sequence files and interface with io.c

QUALITY.CONSTANTS <- list(phred=list(offset=33, min=0, max=93),
                          solexa=list(offset=64, min=-5, max=62),
                          illumina=list(offset=64, min=0, max=62))

NUCLEOTIDES <- c('A', 'T', 'C', 'G', 'N')
NUCLEOTIDES.COLORS <- c('A'='dark green', 'T'='red',
                        'C'='blue', 'G'='black',
                        'N'='purple')

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
  for (i in 1:ncol(qual.dist)) {
    vals <- as.integer(colnames(qual.dist[-1]))
    means[i] <- weighted.mean(vals, w=qual.dist[i, -1])
  }
  if (is.null(length.dist))
    return(mean(means))
  return(weighted.mean(means, w=lengths2weights(length.dist)))
}

readSeqFile <-
# Use the C function summarize_file to create matrices of base
# and quality counts, per position along all reads.
# TODO: Currently only the illumina quality setting has been tested.
function(filename, type='fastq', max.length=400, quality='illumina', hash=TRUE, verbose=FALSE) {
  if (!file.exists(filename))
    stop(sprintf("file '%s' does not exist", filename))

  if (type == 'fasta') {
    qtype <- -1
    quality <- NULL
  } else
    qtype <- which(names(QUALITY.CONSTANTS) == quality) - 1
  
  out <- .Call('summarize_file', filename,
               as.integer(max.length),
               as.integer(qtype),
               as.logical(hash),
               as.logical(verbose))
  
  names(out) <- c('base.freqs', 'seq.lengths', 'qual.freqs')

  obj <- new("SequenceSummary")
  obj@seq.lengths <- .trimArray(out$seq.lengths[-1]) # from C call, this is 0-indexed
  
  if (hash) {
    names(out)[4] <- 'hash'
    obj@hash <- sortSequenceHash(out$hash)
  }

  ## Data cleaning
  obj@base.freqs <- local({
    tmp <- as.data.frame(t(.trimRightCols(out$base.freqs)))
    tmp <- cbind(1:nrow(tmp), tmp)
    colnames(tmp) <- c('position', NUCLEOTIDES)
    return(tmp)})

  if (type == 'fastq') {
    obj@qual.freqs <- local({
      tmp <- .trimRightCols(out$qual.freqs)
      tmp <- as.data.frame(t(.setQualityNames(tmp, quality)))
      tmp <- cbind(1:nrow(tmp), tmp)
      colnames(tmp)[1] <- 'position'
      return(tmp)})
    obj@mean.qual <- meanFromBins(obj@qual.freqs, obj@seq.lengths)
  }
  
  obj@base.props <- local({
    tmp <- prop.table(as.matrix(obj@base.freqs[, -1]), margin=1)
    tmp <- melt(tmp, id='position')
    colnames(tmp) <- c('position', 'base', 'proportion')
    return(tmp)
  })

  obj@hashed <- hash
  obj@quality <- quality
  obj@type <- type
  return(obj)
}
