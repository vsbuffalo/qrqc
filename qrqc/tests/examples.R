## examples.R - examples to ensure everything is functioning correctly

dyn.load('../src/io.so')

source('../R/scan.R')

s <- summarizeFastq('../data/test.fastq', hash=TRUE, verbose=TRUE)

plotQuals(s)
plotBaseProps(s)
plotBaseFreqs(s)
