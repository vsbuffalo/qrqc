## show-methods.R - show methods used in qrqc

## Methods for FASTQSummary
setMethod(show, "FASTQSummary",
          function(object) {
            cat(sprintf("Quality Information for: %s\n", object@filename))
            cat(sprintf(" %i sequences", sum(object@seq.lengths)))
            if (object@hashed)
              cat(sprintf(", %i unique", length(object@hash)))
            cat("\n")
            cat(sprintf(" mean quality: %f\n", round(object@mean.qual, 4)))
            cat(sprintf(" min sequence length: %i\n", min(which(object@seq.lengths > 0))))
            cat(sprintf(" max sequence length: %i\n", length(object@seq.lengths))) # assumes .trimArray called
          })

## Methods for FASTASummary
setMethod(show, "FASTASummary",
          function(object) {
            cat(sprintf("Quality Information for: %s\n", object@filename))
            cat(sprintf(" %i sequences", sum(object@seq.lengths)))
            if (object@hashed)
              cat(sprintf(", %i unique", length(object@hash)))
            cat("\n")
            cat(sprintf(" min sequence length: %i\n", min(which(object@seq.lengths > 0))))
            cat(sprintf(" max sequence length: %i\n", length(object@seq.lengths))) # assumes .trimArray called
          })
