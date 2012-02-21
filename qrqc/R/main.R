## main.R - functions to handle sequence files and interface with io.c

QUALITY.CONSTANTS <- list(phred=list(offset=0, min=4, max=60),
                          sanger=list(offset=33, min=0, max=93),
                          solexa=list(offset=64, min=-5, max=62),
                          illumina=list(offset=64, min=0, max=62))

## for new ggplot2 functions
DNA_BASES_N <- c("A", "T", "G", "C", "N")

###### deprecated start -- to remove when old plotting functions are removed ######
NUCLEOTIDES <- c('A', 'T', 'C', 'G', 'N', 'R', 'Y', 'S',
                 'W', 'K', 'M', 'B', 'D', 'H', 'V', '-')
NUCLEOTIDES.COLORS <- c('A'='dark green', 'T'='red',
                        'C'='blue', 'G'='black',
                        'N'='purple')
other.iupac <- c('R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', '-')


# These are from RColorBrewer, but I include them manually to remove a
# depedency for one call.
other.iupac.colors <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", 
                        "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5")
names(other.iupac.colors) <- other.iupac
NUCLEOTIDES.COLORS <- c(NUCLEOTIDES.COLORS, other.iupac.colors)
###### deprecated end ######

readSeqFile <-
# Use the C function summarize_file to create matrices of base
# and quality counts, per position along all reads.
function(filename, type='fastq', max.length=1000, quality='illumina', hash=TRUE, hash.prop=0.1, verbose=FALSE) {
  if (!file.exists(filename))
    stop(sprintf("file '%s' does not exist", filename))

  if (type == 'fasta') {
    qtype <- -1
    quality <- NULL
  } else {
    tmp <- strsplit(filename, '.', fixed=TRUE)[[1]]
    if (tmp[length(tmp)] == 'fasta') {
      ans <- readline(sprintf("%s's extension is '.fasta'.\nAre you sure you want type='fastq'? [y/n] ", filename))
      if (ans == "n")
        return()
    }
    qtype <- which(names(QUALITY.CONSTANTS) == quality) - 1
  }
  
  out <- .Call('summarize_file', filename,
               as.integer(max.length),
               as.integer(qtype),
               as.logical(hash),
               as.numeric(hash.prop),
               as.logical(verbose))
  
  names(out) <- c('base.freqs', 'seq.lengths', 'qual.freqs')

  class.types <- c('fasta'="FASTASummary", 'fastq'="FASTQSummary")
  obj <- new(class.types[type], filename=filename)
  
  obj@seq.lengths <- .trimArray(out$seq.lengths[-1]) # from C call, this is 0-indexed
  
  if (hash) {
    names(out)[4] <- 'hash'
    obj@hash <- sortSequenceHash(out$hash)
    obj@hash.prop <- hash.prop
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
    obj@quality <- quality
  }
  obj@hashed <- hash
  return(obj)
}
