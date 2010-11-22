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
    stop('Conversions to Sanger qualities only supported')

  
  
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
.Call('convert_quality', "A")
