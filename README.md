# HHGA

![orb](https://raw.githubusercontent.com/ekg/hhga/master/images/orb.jpg)

## Have Haplotypes Genotypes Alleles

### A tool to prepare genome alignments and genotypes for consumption by the Vowpal Wabbit

[![Build Status](https://travis-ci.org/ekg/hhga.svg)](https://travis-ci.org/ekg/hhga)

## Overview

HHGA transforms genomic alignments and candidate haplotypes into a feature-space representation suitable for use by machine learning systems, in particular [Vowpal Wabbit](https://github.com/JohnLangford/vowpal_wabbit).
We take inspiration from genome visualizers that are often used to inspect candidate variant calls.
For example, `samtools tview` implements a text-based view of alignments in a particular region of the reference genome with encodings representing the read placement, direction, non-reference bases in the alignments, gaps (insertions and deletions are both shown faithfully using a gap character), and the qualities of mismatches. 


![tview](https://raw.githubusercontent.com/ekg/hhga/master/images/tview.png)

Additional tools provide mechanisms to convert these representations into human-readable formats (both text and image) for debugging.

HHGA can be thought of as an example decision synthesizer.
We use it in combination with genomic truth sets to generate examples from sequencing data that can be used in decision learning.
Our primary interest is to approximate P(genotype|alignments), where the genotype is composed of N haplotypes (typically 2), and the alignments represent the observational data we have at the genomic locus.

## Feature space representation

We represent the MSA of the reads and the reference, the estimated per-base errors, the read strands releative to the reference, the alignment mapping qualities, and the candidate haplotypes in the genotype of the sample at the locus using the namespaced libSVM format used by Vowpal Wabbit.

As an example, if this is a tview-like representation of a set of alignments at a putative SNP site, augmented with candidate haplotypes. Generated via `hhga -b minigiab/9251-9252.bam -f minigiab/q.fa -v minigiab/h.vcf.gz -r q:9251-9252 -w 20 -t` from the `test/` directory.

```txt
reference   ATTGTGCCAAGTTCTTTCTT
hap                   GTTCT     
hap                   G----     
s    xypzi  AT                   60 chr22.bin8.cram:166:8383
s    xypzi  ATTGT                60 chr22.bin8.cram:166:8410
 o   xypzi  ATTGTGCCAAGTTCTTTC   60 chr22.bin8.cram:166:8436
 o  f ypzi  ATTGTGCCAAG----TTCTT 60 chr22.bin8.cram:166:8455
s   f ypzi  ATTGTGCCAAG----TTCTT 60 chr22.bin8.cram:166:8390
 o  f ypzi  ATTGTGCCAAGTTCTTTCTT 60 chr22.bin8.cram:166:8457
s   f ypzi  ATTGTGCCAAGTTCTTTCTT 60 chr22.bin8.cram:166:8436
s    xypzi  ATTGTGCCAAG----TTCTT 60 chr22.bin8.cram:166:8400
 o   xypzi  ATTGTGCCAAG----TTCTT 60 chr22.bin8.cram:166:8460
 o   xypzi       GCCAAGTTCTTTCTT 60 chr22.bin8.cram:166:8475
s   f ypzi                   CTT 60 chr22.bin8.cram:166:8395
```

We use six symbols to encode the MSA, `{ A, T, G, C, N, U, M }`, where `N` is the degenerate symbol (it could represent any of A, T, G, or C), `U` is a gap symbol required to normalize the MSA into a matrix, and `M` is a symbol indicating if we don't have any information at the position in the MSA.

Various features of the reads are represented in other namespaces.

The feature space transformation of this site would be given by `hhga -b minigiab/9251-9252.bam -f minigiab/q.fa -v minigiab/h.vcf.gz -r q:9251-9252 -w 20 -c 1`. The output will come on a single line.

```txt
1 |ref A:1 T:1 T:1 G:1 T:1 G:1 C:1 C:1 A:1 A:1 G:1 T:1 T:1 C:1 T:1 T:1 T:1 C:1 T:1 T:1 |hap0 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 G:1 T:1 T:1 C:1 T:1 M:1 M:1 M:1 M:1 M:1 |hap1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 G:1 U:1 U:1 U:1 U:1 M:1 M:1 M:1 M:1 M:1 |aln0 A:0.9998 T:0.999499 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 |aln1 A:0.9998 T:0.9998 T:0.9998 G:0.999499 T:0.999499 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 |aln2 A:0.9999 T:0.9999 T:0.99369 G:0.9998 T:0.9998 G:0.9998 C:0.9998 C:0.9998 A:0.9999 A:0.9998 G:0.9998 T:0.998005 T:0.999499 C:0.9999 T:0.9999 T:0.9998 T:0.9999 C:0.9998 M:1 M:1 |aln3 A:0.9999 T:0.9999 T:0.9999 G:0.9999 T:0.9999 G:0.9999 C:0.9999 C:0.9999 A:0.9999 A:0.999499 G:0.9999 U:0.9999 U:0.999499 U:0.9999 U:0.9999 T:0.9999 T:0.9999 C:0.9999 T:0.9999 T:0.9999 |aln4 A:0.9999 T:0.998005 T:0.9998 G:0.9999 T:0.9999 G:0.9999 C:0.9999 C:0.9999 A:0.9998 A:0.9999 G:0.9999 U:0.9998 U:0.9999 U:0.9999 U:0.9999 T:0.9999 T:0.9999 C:0.9999 T:0.9999 T:0.9999 |aln5 A:0.9999 T:0.9999 T:0.999499 G:0.999499 T:0.9999 G:0.9999 C:0.9999 C:0.9999 A:0.968377 A:0.999499 G:0.998005 T:0.999499 T:0.968377 C:0.9998 T:0.9999 T:0.9999 T:0.999499 C:0.9999 T:0.99369 T:0.9999 |aln6 A:0.9998 T:0.999499 T:0.9999 G:0.9998 T:0.9999 G:0.9999 C:0.9999 C:0.9999 A:0.9999 A:0.9999 G:0.9999 T:0.9999 T:0.9999 C:0.9998 T:0.9999 T:0.9999 T:0.9999 C:0.9999 T:0.9999 T:0.9999 |aln7 A:0.9999 T:0.9999 T:0.9998 G:0.9998 T:0.9999 G:0.9998 C:0.9999 C:0.9998 A:0.9998 A:0.9999 G:0.9998 U:0.9998 U:0.9999 U:0.9998 U:0.9999 T:0.9999 T:0.9998 C:0.9999 T:0.9998 T:0.9999 |aln8 A:0.9999 T:0.9999 T:0.9999 G:0.9999 T:0.9999 G:0.9999 C:0.9999 C:0.9999 A:0.9999 A:0.9999 G:0.9999 U:0.9999 U:0.9999 U:0.9999 U:0.9999 T:0.9999 T:0.9999 C:0.9999 T:0.9999 T:0.9999 |aln9 M:1 M:1 M:1 M:1 M:1 G:0.999499 C:0.9998 C:0.9998 A:0.9998 A:0.9998 G:0.9998 T:0.9999 T:0.9999 C:0.9999 T:0.9999 T:0.9999 T:0.9999 C:0.9999 T:0.9999 T:0.9999 |aln10 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 C:0.9998 T:0.9998 T:0.9998 |mapq aln0:0.999999 aln1:0.999999 aln10:0.999999 aln2:0.999999 aln3:0.999999 aln4:0.999999 aln5:0.999999 aln6:0.999999 aln7:0.999999 aln8:0.999999 aln9:0.999999 |strand aln0:1 aln1:1 aln10:1 aln4:1 aln6:1 aln7:1 |ostrand aln2:1 aln3:1 aln5:1 aln8:1 aln9:1 |dup |qcfail |fmate aln10:1 aln3:1 aln4:1 aln5:1 aln6:1 |xmate aln0:1 aln1:1 aln2:1 aln7:1 aln8:1 aln9:1 |ymap aln0:1 aln1:1 aln10:1 aln2:1 aln3:1 aln4:1 aln5:1 aln6:1 aln7:1 aln8:1 aln9:1 |paired aln0:1 aln1:1 aln10:1 aln2:1 aln3:1 aln4:1 aln5:1 aln6:1 aln7:1 aln8:1 aln9:1 |zprimary aln0:1 aln1:1 aln10:1 aln2:1 aln3:1 aln4:1 aln5:1 aln6:1 aln7:1 aln8:1 aln9:1 |iproper aln0:1 aln1:1 aln10:1 aln2:1 aln3:1 aln4:1 aln5:1 aln6:1 aln7:1 aln8:1 aln9:1 
```

For viewing we can pipe the hhga format output through `sed s/\|/\\n\|/g` to convert the spaces into newlines.

```txt
1 
|ref A:1 T:1 T:1 G:1 T:1 G:1 C:1 C:1 A:1 A:1 G:1 T:1 T:1 C:1 T:1 T:1 T:1 C:1 T:1 T:1 
|hap0 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 G:1 T:1 T:1 C:1 T:1 M:1 M:1 M:1 M:1 M:1 
|hap1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 G:1 U:1 U:1 U:1 U:1 M:1 M:1 M:1 M:1 M:1 
|aln0 A:0.9998 T:0.999499 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 
|aln1 A:0.9998 T:0.9998 T:0.9998 G:0.999499 T:0.999499 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 
|aln2 A:0.9999 T:0.9999 T:0.99369 G:0.9998 T:0.9998 G:0.9998 C:0.9998 C:0.9998 A:0.9999 A:0.9998 G:0.9998 T:0.998005 T:0.999499 C:0.9999 T:0.9999 T:0.9998 T:0.9999 C:0.9998 M:1 M:1 
|aln3 A:0.9999 T:0.9999 T:0.9999 G:0.9999 T:0.9999 G:0.9999 C:0.9999 C:0.9999 A:0.9999 A:0.999499 G:0.9999 U:0.9999 U:0.999499 U:0.9999 U:0.9999 T:0.9999 T:0.9999 C:0.9999 T:0.9999 T:0.9999 
|aln4 A:0.9999 T:0.998005 T:0.9998 G:0.9999 T:0.9999 G:0.9999 C:0.9999 C:0.9999 A:0.9998 A:0.9999 G:0.9999 U:0.9998 U:0.9999 U:0.9999 U:0.9999 T:0.9999 T:0.9999 C:0.9999 T:0.9999 T:0.9999 
|aln5 A:0.9999 T:0.9999 T:0.999499 G:0.999499 T:0.9999 G:0.9999 C:0.9999 C:0.9999 A:0.968377 A:0.999499 G:0.998005 T:0.999499 T:0.968377 C:0.9998 T:0.9999 T:0.9999 T:0.999499 C:0.9999 T:0.99369 T:0.9999 
|aln6 A:0.9998 T:0.999499 T:0.9999 G:0.9998 T:0.9999 G:0.9999 C:0.9999 C:0.9999 A:0.9999 A:0.9999 G:0.9999 T:0.9999 T:0.9999 C:0.9998 T:0.9999 T:0.9999 T:0.9999 C:0.9999 T:0.9999 T:0.9999 
|aln7 A:0.9999 T:0.9999 T:0.9998 G:0.9998 T:0.9999 G:0.9998 C:0.9999 C:0.9998 A:0.9998 A:0.9999 G:0.9998 U:0.9998 U:0.9999 U:0.9998 U:0.9999 T:0.9999 T:0.9998 C:0.9999 T:0.9998 T:0.9999 
|aln8 A:0.9999 T:0.9999 T:0.9999 G:0.9999 T:0.9999 G:0.9999 C:0.9999 C:0.9999 A:0.9999 A:0.9999 G:0.9999 U:0.9999 U:0.9999 U:0.9999 U:0.9999 T:0.9999 T:0.9999 C:0.9999 T:0.9999 T:0.9999 
|aln9 M:1 M:1 M:1 M:1 M:1 G:0.999499 C:0.9998 C:0.9998 A:0.9998 A:0.9998 G:0.9998 T:0.9999 T:0.9999 C:0.9999 T:0.9999 T:0.9999 T:0.9999 C:0.9999 T:0.9999 T:0.9999 
|aln10 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 M:1 C:0.9998 T:0.9998 T:0.9998 
|mapq aln0:0.999999 aln1:0.999999 aln10:0.999999 aln2:0.999999 aln3:0.999999 aln4:0.999999 aln5:0.999999 aln6:0.999999 aln7:0.999999 aln8:0.999999 aln9:0.999999 
|strand aln0:1 aln1:1 aln10:1 aln4:1 aln6:1 aln7:1 
|ostrand aln2:1 aln3:1 aln5:1 aln8:1 aln9:1 
|dup 
|qcfail 
|fmate aln10:1 aln3:1 aln4:1 aln5:1 aln6:1 
|xmate aln0:1 aln1:1 aln2:1 aln7:1 aln8:1 aln9:1 
|ymap aln0:1 aln1:1 aln10:1 aln2:1 aln3:1 aln4:1 aln5:1 aln6:1 aln7:1 aln8:1 aln9:1 
|paired aln0:1 aln1:1 aln10:1 aln2:1 aln3:1 aln4:1 aln5:1 aln6:1 aln7:1 aln8:1 aln9:1 
|zprimary aln0:1 aln1:1 aln10:1 aln2:1 aln3:1 aln4:1 aln5:1 aln6:1 aln7:1 aln8:1 aln9:1 
|iproper aln0:1 aln1:1 aln10:1 aln2:1 aln3:1 aln4:1 aln5:1 aln6:1 aln7:1 aln8:1 aln9:1 
```

The first entry in the line defines the class of the example. By convention, we say the example is 1 if the haplotypes are correct, and -1 otherwise.

* `ref` : the reference
* `hap0`: the first haplotype
* `hap1`: the second haplotype
* `aln*` : the alignments
* `mapq` : p(mapping correct)
* `strand` : 1 if reversed, 0 otherwise (not written)
* `ostrand` : mate's strand, 1 if reversed
* `dup` : 1 if the read is marked as a duplicate
* `qcfail` : 1 if the read has failed some QC
* `fmate` : 1 if the read is the first mate
* `xmate` : 1 if the read is the second mate
* `ymap` : 1 if the mate is mapped (this mate must be mapped to be output)
* `paired` : 1 if the read is paired
* `zprimary` : 1 if the read is the "primary" alignment
* `iproper` : 1 if the read is in a "proper" pair, in the expected configuration and distance from its mate


## Usage

Sketch of usage. We extract the feature space representation using windows of size 100 at sites that are defined in the candidates in the variant input file. We use `truth.vcf.gz` to indicate which genotypes and alleles are true. All others in `vars.vcf.gz` are assumed false.

```bash
( hhga -v trues.vcf.gz  -w 100 -b aln.bam -f ref.fa -c 1 \
  hhga -v falses.vcf.gz -w 100 -b aln.bam -f ref.fa -c -1 ) \
    | shuf \
    | vw --save_resume  --ngram 5 --skips 3 --loss_function logistic --interactions rhmsa
```

## Building

```
git clone --recursive https://github.com/ekg/hhga.git
cd hhga
source ./source_me.sh
make
make test
```

The `hhga` executable is at `bin/hhga`.

## TODO

Make this go, add a method to generate more training data per site by varying the length of the candidate haplotypes (e.g. full-window vs. just SNP view), add image output for debugging.
