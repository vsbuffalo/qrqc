## AllClasses.R - classes used in qrqc

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

