# HHGA

![orb](https://raw.githubusercontent.com/ekg/hhga/master/images/orb.jpg)

## Have Haplotypes Genotypes Alleles

### A tool to prepare genome alignments and genotypes for consumption by the Vowpal Wabbit

## Overview

HHGA transforms genomic alignments and candidate haplotypes into a feature-space representation suitable for use by machine learning systems. We take inspiration from genome visualizers that are often used to inspect candidate variant calls. For example, `samtools tview` implements a text-based view of alignments in a particular region of the reference genome with encodings representing the read placement, direction, non-reference bases in the alignments, gaps (insertions and deletions are both shown faithfully using a gap character), and the qualities of mismatches. Note that the tview representation is matrix-like, but the default behavior of HHGA does not generate a true matrix.

![tview](https://raw.githubusercontent.com/ekg/hhga/master/images/tview.png)

Additional tools provide mechanisms to convert these representations into human-readable formats (both text and image) for debugging.

HHGA can be thought of as an example decision synthesizer. We use it in combination with genomic truth sets to generate examples from sequencing data that can be used in decision learning. Our primary interest is to approximate P(genotype|alignments), where the genotype is composed of N haplotypes (typically 2), and the alignments represent the observational data we have at the genomic locus.

## Feature space representation

We represent the MSA of the reads and the reference, the estimated per-base errors, the read strands releative to the reference, the alignment mapping qualities, and the candidate haplotypes in the genotype of the sample at the locus using the namespaced libSVM format used by Vowpal Wabbit.

As an example, if this is a tview-like representation of a set of alignments at a putative SNP site, augmented with candidate haplotypes:

```txt
ref:  TAT
hap1: ...
hap2: .G.
aln1: ,,,
aln2: .G.
aln3:  g,
aln4: ..
```

Assume our mapping qualities for the alignments are aln1: 60, aln2: 20, aln3: 50, aln4: 30.

We use six symbols to encode the MSA, `{ A, T, G, C, N, U, M }`, where `N` is the degenerate symbol (it could represent any of A, T, G, or C), `U` is a gap symbol required to normalize the MSA into a matrix, and `M` is a symbol indicating if we don't have any information at the position in the MSA.

Assuming the genotype is true and there are a mixture of 3 different error probabilities for the bases in the reads, the feature space transformation of this site would be:

```txt
1 |ref    1T:1 2A:1 3T:1
  |hap1   1T:1 2A:1 3T:1
  |hap2   1T:1 2G:1 3T:1
  |mapq   aln1:60 aln2:20 aln3:50 aln4:30
  |strand aln1:1  aln2:0  aln3:1  aln4:0
  |aln1   1A:0.02 1T:0.94 1G:0.02 1C:0.02 2A:0.91 2T:0.03 2G:0.03 2C:0.03 3A:0.003 3T:0.991 3G:0.003 3C:0.003
  |aln2   1A:0.02 1T:0.94 1G:0.02 1C:0.02 2A:0.003 2T:0.003 2G:0.991 2C:0.003 3A:0.02 3T:0.94 3G:0.02 3C:0.02
  |aln3   1M:1 2A:0.02 2T:0.02 2G:0.94 2C:0.02 3A:0.02 3T:0.94 3G:0.02 3C:0.02
  |aln4   1A:0.02 1T:0.94 1G:0.02 1C:0.02 2A:0.02 2T:0.02 2G:0.94 2C:0.02 3M:1
```

Newlines are included here only for legibility. This would be on one line of the output of HHGA.

As in libSVM format, the first entry in the output defines the class of the example. By convention, we say the example is 1 if the haplotypes are correct, and 0 otherwise.

## Usage

Sketch of usage. We extract the feature space representation using windows of size 100 at sites that are defined in the candidates in the variant input file. We use `truth.vcf.gz` to indicate which genotypes and alleles are true. All others in `vars.vcf.gz` are assumed false.

```bash
hhga -t truth.vcf.gz -v vars.vcf.gz -w 100 -b aln.bam -f ref.fa >examples.hhga
```

The output can now be used in [Vowpal Wabbit](https://github.com/JohnLangford/vowpal_wabbit) and other systems which support the same input format.

## TODO

Make this go, add a method to generate more training data per site by varying the length of the candidate haplotypes (e.g. full-window vs. just SNP view), add image output for debugging.
