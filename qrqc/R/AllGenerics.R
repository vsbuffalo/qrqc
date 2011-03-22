## AllGenerics.R 

setGeneric("plotQuals", signature="obj", function(obj, ylim='relative', lowess=TRUE, histogram=TRUE, legend=TRUE) standardGeneric("plotQuals"))
setGeneric("plotBases", signature="obj", function(obj, type="freq", bases=NULL, legend=TRUE) standardGeneric("plotBases"))
setGeneric("plotSeqLengths", signature="obj", function(obj) standardGeneric("plotSeqLengths"))
setGeneric("plotGC", signature="obj", function(obj) standardGeneric("plotGC"))
setGeneric("makeReport", signature="obj", function(obj, outputDir=NULL) standardGeneric("makeReport"))
