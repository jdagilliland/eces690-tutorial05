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
IDBA performs better than similar earlier tools (e.g. SSAKE, VCAKE,
SHARCGS) in terms of memory use and time.

### Using IDBA

For assembling genomic data without a reference, the IDBA package recommends
using IDBA-UD (installed on Proteus as `idba_ud`).

## Read Mapping

The most common tool for read mapping is bowtie2.

### Using bowtie2

#### Building an index

When using bowtie2, it is crucial to have an index for the reference
genome/transcriptome.
This reference genome can come from a pre-built index, or can be built from
a set of FASTA files that detail the genes in question.
To build an index, bowtie2 provides the `bowtie2-build` command, which
takes as arguments either a single FASTA file or a comma separated list
thereof.
The index allows the read mapper to efficiently decide the location to which
to map a given read.
To build an index, however, already requires that you have FASTA file(s)
detailing the genes to which you want the reads to map.

The following is an example of how such a reference index would be built:
`bowtie2-build reference_genome.fasta refgen`

In this case, `refgen` is a name that we choose arbitrarily that tags all of
the files that constitute the index.
`reference_genome.fasta` must point to the actual FASTA file that contains
our reference genome.

#### Performing the alignment

Performing the alignment requires the actual reads (in fastq format), as
well as the index that was built in the last step.
It outputs a SAM file, that details the alignment results.

Consider, for example, the following:
`bowtie2 -x refgen -U reads_1.fastq -S eg1.sam`


