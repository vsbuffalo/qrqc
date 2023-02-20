# qrqc - Quick Read Quality Control

qrqc and all supporting documentation 
Copyright (c) Vince Buffalo, 2011-2012

Contact: Vince Buffalo <vsbuffaloAAAAA@gmail.com> (with the poly-A tail removed)

If you wish to report a bug, please open an issue on Github
(http://github.com/vsbuffalo/qrqc/issues) or post it on the
Bioconductor support site (https://support.bioconductor.org/).
You can contact me personally as well, but please open an issue first.

## About

qrqc (short for "Quick Read Quality Control") is a fast and extensible
package that reports basic quality and summary statistics on FASTQ and
FASTA files, including base and quality distribution by position,
sequence length distribution, and common sequences.

## License

GNU General Public License, version 2.

## FAQ

### Why `ggplot2`?

I've had some feature requests for `qrqc` since its release, mostly
related to customizing the graphics. Since data accessibility and
custom graphics were the reason I created `qrqc`, I initially rewrote
`qrqc` to provide more graphics options through `lattice`. However,
all the graphics parameters I added led to large numbers of arguments
to functions and high complexity. This rewrite uses `ggplot2`, which
is a very excellent way to create graphics as any graphics object can
be further manipulated.

### Why do you use Monte Carlo simulations to generate the smooth curve?

`qrqc` is fast because it bins the quality scores of bases by
positions; there is data summarization done by `readSeqFile`. To
create a smooth curve, the function needs multiple data points (not
binned data), which I simulate via Monte Carlo draws from the quality
distribution by position. This is an approximation, but it leads to a
smooth curve which can create a useful visual tool in assessing
quality drops.

### What do I do about bad quality regions?

Illumina reads often have poor 3'-end qualities. I've noticed that
HiSeq machines also produce poor quality 5'-ends. For increased
mapping rates and better assemblies, it is generally advisable that
these poor quality regions be trimmed off. Nik Joshi's took `sickle`
tool can do this; you can get it here
<http://github.com/najoshi/sickle>.

3'-end adapter contamination can be difficult to recognize (and thus
remove) due to poor quality and likely incorrect bases. I've developed
a tool called `scythe` that removes 
