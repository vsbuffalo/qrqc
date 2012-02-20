## AllGenerics.R 

setGeneric("plotQuals", signature="obj", function(obj, ylim='relative', lowess=TRUE, histogram=TRUE, legend=TRUE) standardGeneric("plotQuals"))
setGeneric("plotBases", signature="obj", function(obj, type="freq", bases=NULL, legend=TRUE) standardGeneric("plotBases"))
setGeneric("plotSeqLengths", signature="obj", function(obj) standardGeneric("plotSeqLengths"))
setGeneric("plotGC", signature="obj", function(obj) standardGeneric("plotGC"))
setGeneric("makeReport", signature="obj", function(obj, outputDir=".") standardGeneric("makeReport"))


## ggplot2 generics
setGeneric("qualPlot", signature="x",
           function(x, smooth=TRUE, extreme.color="grey", quartile.color="orange",
                    mean.color="blue", median.color=NULL) standardGeneric("qualPlot"))
setGeneric("gcPlot", signature="x", function(x, color="red") standardGeneric("gcPlot"))
setGeneric("basePlot", signature="x", function(x, geom=c("line", "bar", "dodge"), type=c("frequency", "proportion"),
                         bases=DNA_BASES_N, colorvalues=getBioColor("DNA_BASES_N")) standardGeneric("basePlot"))
setGeneric("seqlenPlot", signature="x", function(x) standardGeneric("seqlenPlot"))
setGeneric("getQual", signature="x", function(x) standardGeneric("getQual"))
setGeneric("getGC", signature="x", function(x) standardGeneric("getGC"))
setGeneric("getBase", signature="x", function(x, drop=TRUE) standardGeneric("getBase"))
setGeneric("getBaseProp", signature="x", function(x, drop=TRUE) standardGeneric("getBaseProp"))
setGeneric("getSeqlen", signature="x", function(x) standardGeneric("getSeqlen"))
setGeneric("getMCQual", signature="x", function(x, n=100) standardGeneric("getMCQual"))
