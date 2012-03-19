## AllClasses.R - classes used in qrqc

setClass("SequenceSummary",
         representation=representation(
           "VIRTUAL",
           filename='character',
           base.freqs='data.frame',
           seq.lengths='integer',
           hash='numeric',
           hash.prop='numeric',
           kmer='data.frame',
           k='integer',
           hashed='logical',
           kmers.hashed='logical'))

setClass("FASTASummary",
         contains="SequenceSummary")

setClass("FASTQSummary",
         representation=representation(
           qual.freqs='data.frame',
           mean.qual='numeric',
           quality='character'),
         contains="SequenceSummary")

