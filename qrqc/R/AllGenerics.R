## AllGenerics.R 

setGeneric("plotQuals", signature="obj", function(obj, ylim='relative', lowess=TRUE, histogram=TRUE) standardGeneric("plotQuals"))
setGeneric("plotBaseProps", signature="obj", function(obj) standardGeneric("plotBaseProps"))
setGeneric("plotBaseFreqs", signature="obj", function(obj, bases=NULL) standardGeneric("plotBaseFreqs"))
setGeneric("plotSeqLengths", signature="obj", function(obj) standardGeneric("plotSeqLengths"))
setGeneric("plotGC", signature="obj", function(obj) standardGeneric("plotGC"))
setGeneric("makeReport", signature="obj", function(obj, filename=NULL) standardGeneric("makeReport"))
