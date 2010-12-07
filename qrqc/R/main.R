## main.R - functions to handle sequence files and interface with io.c

QUALITY.CONSTANTS <- list(phred=list(offset=33, min=0, max=93),
                          solexa=list(offset=64, min=-5, max=62),
                          illumina=list(offset=64, min=0, max=62))

NUCLEOTIDES <- c('A', 'T', 'C', 'G', 'N')
NUCLEOTIDES.COLORS <- c('A'='dark green', 'T'='red',
                        'C'='blue', 'G'='black',
                        'N'='purple')


setClass("SequenceSummary",
         representation=representation(
           base.freqs='data.frame',
           base.props='data.frame',
           qual.freqs='data.frame',
           quality='character',
           seq.lengths='integer',
           type='character'))

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

  if (type != -1) {
    obj@qual.freqs <- local({
      tmp <- .trimRightCols(out$qual.freqs)
      tmp <- as.data.frame(t(.setQualityNames(tmp, quality)))
      tmp <- cbind(1:nrow(tmp), tmp)
      colnames(tmp)[1] <- 'position'
      return(tmp)})
  }
  
  obj@base.props <- local({
    tmp <- prop.table(as.matrix(out$base.freqs[, -1]), margin=1)
    tmp <- melt(tmp, id='position')
    colnames(tmp) <- c('position', 'base', 'proportion')
    return(tmp)
  })

  obj@quality <- quality
  obj@type <- type
  return(obj)
}
