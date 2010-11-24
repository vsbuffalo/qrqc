## scan.R - Scan sequenece files

QUALITY.CONSTANTS <- list(phred=list(offset=33, min=0, max=93),
                          solexa=list(offset=64, min=-5, max=62),
                          illumina=list(offset=64, min=0, max=62))
  
setQualityNames <-
#
function(matrix, quality) {
  constants <- QUALITY.CONSTANTS[[quality]]
  rownames(matrix) <- constants$min:constants$max
  return(matrix)
}

trimRightCols <-
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
function(filename, max.length=100, quality='illumina') {
  if (!file.exists(filename))
    stop(sprintf("file '%s' does not exist", filename))
  hash <- new.env(hash=TRUE)
  out <- .Call('summarize_fastq_file', filename,
               as.integer(max.length),
               hash,
               as.integer(which(names(QUALITY.CONSTANTS) == quality) - 1))

  names(out) <- c('base.freqs', 'qual.freqs')
  out$qual.freqs <- setQualityNames(trimRightCols(out$qual.freqs), quality)
  out$hash <- sapply(ls(envir=hash), function(x) get(x, envir=hash))
  return(out)
}

s <- summarizeFastq('test.fastq')
