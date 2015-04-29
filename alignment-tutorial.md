# Aligning sequencing data

Sequencing data is collected to study the content of either genomes or
transcriptomes.
Generally the way the data presents is as short (<200 bases) reads of
nucleotide sequences.
In order for data so collected to be meaningful, we must learn to which
genes or other genomic or transcriptomic entities they belong.
Thus we can form a quantitative picture of what is going on in the sample.

There are fundamentally two approaches to this.
In *de novo* alignment we treat the collection of nucleic acid sequences
(reads) as the only information we have, and the task set before us is to
determine how these reads fit together and overlap like pieces of a puzzle.
The challenges here include:
* huge number of reads
* sequence homology (ambiguous alignment)
In read mapping, a reference genome or transcriptome is known in addition to
the reads, and we simply map reads to where they line up with the reference data.

## De Novo Alignment

De novo alignment may be carried out using the IDBA package.
IDBA stands for Iterative De Bruijn Assembler.
It uses De Bruijn graphs to align sequences even in the absence of
a reference.

## Read Mapping

The most common tool for read mapping is bowtie2.
