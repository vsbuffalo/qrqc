## oo.R - classes and generic functions

setClass("SequenceSummary",
         representation=representation(
           base.freqs='data.frame',
           base.props='data.frame',
           qual.freqs='data.frame',
           mean.qual='numeric',
           quality='character',
           seq.lengths='integer',
           type='character',
           hash='integer',
           hashed='logical'))

setMethod(show, "SequenceSummary",
          function(object) {
            cat("SequenceSummary object\n")
            cat(sprintf(" %i sequences, ", sum(object@seq.lengths)))
            if (object@hashed)
              cat(sprintf(" %i unique", length(object@hash)))
            cat("\n")
            cat(sprintf(" mean quality: %f\n", round(object@mean.qual, 4)))
            cat(sprintf(" min sequence length: %i\n", min(which(object@seq.lengths > 0))))
            cat(sprintf(" max sequence length: %i\n", length(object@seq.lengths))) # assumes .trimArray called
          })


