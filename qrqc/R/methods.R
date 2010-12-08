## methods.R - methods used in qrqc

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


