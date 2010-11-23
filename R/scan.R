## scan.R - Scan sequenece files

initReadQC <-
# Initialize a read QC storage object. `seq.length` is for
# pre-allocation and is the maximum read/quality length.
function(seq.length=100) {
  quals <- vector('list', seq.length)
  bases <- vector('list', seq.length)

  obj <- list(raw.quals=quals, raw.bases=bases)
  class(obj) <- 'ReadQC'
}

convertQuality <-
#
function(qual, to='sanger') {
  if (tolower(to) != 'sanger')
    stop("Conversions to Sanger qualities only supported")
  if (any(!is.character(qual)))
    stop("qualities must be character strings")
  
  if (length(qual) == 1)
    return(.Call('string_to_base_qualities', qual))
  else
    return(lapply(qual, function(q) .Call('string_to_base_qualities', q)))
}

scanFASTQ <-
#
function(filename, obj) {
  
}

scanFASTA <-
#
function(filename, obj) {

}

scanReads <-
# Read in bases, updating a ReadQC object's raw entries
function(filename, type='fastq') {
  obj <- initReadQC()
  
}

dyn.load('io.so')
.Call('string_to_base_qualities', "ABC")

dyn.load('io.so')
base.counts <- vector('list', 100)
base.counts <- lapply(base.counts, function(x) integer(5))
.Call('summarize_fastq_file', 'test.fastq', base.counts)
