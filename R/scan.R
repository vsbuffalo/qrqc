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

  if (hash)
    names(out)[3] <- 'hash'

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
    tmp <- prop.table(as.matrix(obj$base.freqs[, -1]), margin=1)
    tmp <- melt(tmp, id='position')
    colnames(tmp) <- c('position', 'base', 'proportion')
    return(tmp)
  })
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
