# HPSEQ - high-performance next-gen sequencing tools #

Some fast tools for NGS data.

## Available tools: ##

* hpcut - Simple adapter trimming based on Hamming distance.
* hpscan_cw - like grep for (pairs of) FASTQ files (GNU fgrep actually).
* hpscan_sw - like the above but does full Smith-Waterman alignment

## INSTALLATION ##

HPSEQ tools have build time dependencies on GNU grep, klib, and htslib.
The dependencies can be downloaded and built if required using `make deps`.

```bash
git clone "https://github.com/mcieslik-mctp/hpseq.git"
cd hpseq
make deps
make

```
