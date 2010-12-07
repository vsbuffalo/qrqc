## oo.R - classes and generic functions

setClass("SequenceSummary",
         representation=representation(
           base.freqs='data.frame',
           base.props='data.frame',
           qual.freqs='data.frame',
           quality='character',
           seq.lengths='integer',
           type='character'))

setMethod(show, "SequenceSummary",
          function(object) {
            cat("SequenceSummary object\n")
          })


