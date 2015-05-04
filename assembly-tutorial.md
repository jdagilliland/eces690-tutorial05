# Aligning sequencing data

Sequencing data is collected to study the content of either genomes or
transcriptomes.
Generally the way the data presents is as short (<200 bases) reads of
nucleotide sequences.
In order for data so collected to be meaningful, we must learn to which
genes or other genomic or transcriptomic entities they belong.
Thus we can form a quantitative picture of what is going on in the sample.

There are fundamentally two major steps to this.
In read mapping (the last step), a reference genome or transcriptome is
known in addition to
the reads, and we simply map reads to where they line up with the reference data.
The challenges here include:
* doing it faster than BLAST

Read mapping presupposes the availability of a reference genome.
This is where *de novo* alignment comes into play.
In *de novo* alignment we treat the collection of nucleic acid sequences
(reads) as the only information we have, and the task set before us is to
determine how these reads fit together and overlap like pieces of a puzzle.
The challenges here include:
* huge number of reads
* imperfect read accuracy
* sequence homology (ambiguous alignment)

## *De Novo* Alignment

There are numerous available tools for *de novo* alignment, most of which
use either the string graph algorithm, or the De Bruijn graph algorithm,
both of which construct data structures that link parts of the overall
sequence together by identifying overlapping reads.
The tool we will use (IDBA) is an example of a De Bruijn graph method.
IDBA stands for Iterative De Bruijn Assembler.
It uses De Bruijn graphs to align sequences even in the absence of
a reference.
It is iterative in the sense that it uses different length $k$-mers to
construct the graph.
IDBA performs better than similar earlier tools (e.g. SSAKE, VCAKE,
SHARCGS) in terms of memory use and time.

### Using IDBA

For assembling genomic data without a reference, the IDBA package recommends
using IDBA-UD (installed on Proteus as `idba_ud`).
Since our goal here is to create a reference relevant to a similar set of
samples to which we will later align the reads of each sample, we use all of
the samples we hope to analyze as sources for the assembly, in order to
improve read depth.

For this, something along the lines of:
`idba_ud -l *.fastq -o outputdir`
is exactly what we need.
It uses all of the reads from each fastq file, and determines contigs for
each $k$ value used, placing them in `outputdir`.
The `-l` flag is used to signal that we are providing the algorithm with
reads longer than 128 bases (in our case: 200 bases).

A number of other options exist for tuning the algorithm.
For greater detail see [the help output](idba_ud.help), which can also be
obtained by `idba_ud --help`.

Run with default parameters (as we did), `idba_ud` outputs a number of
files.

## Read Mapping

One of the most common tools for read mapping is bowtie2.

### Using bowtie2

#### Building an index

When using bowtie2, it is crucial to have an index for the reference
genome/transcriptome.
This reference genome can come from a pre-built index, or can be built from
a set of FASTA files that contain the contiguous sequences to which reads
are meant to be aligned.
To build an index, bowtie2 provides the `bowtie2-build` command, which
takes as arguments either a single FASTA file or a comma separated list
thereof.
The index allows the read mapper to efficiently decide the location to which
to map a given read.
To build an index, however, already requires that you have FASTA file(s)
detailing the genes to which you want the reads to map.

The following is an example of how such a reference index would be built:
`bowtie2-build reference_genome.fa,and_another.fa refgen`

In this case, `refgen` is a name that we choose arbitrarily that tags all of
the files that constitute the index.
`reference_genome.fa` must point to the actual FASTA file that contains
part of our reference genome.
Additional files containing parts of our reference genome must be provided
separated by commas.

#### Performing the alignment

Performing the alignment requires the actual reads (in fastq format), as
well as the index that was built in the last step.
It outputs a SAM file, that details the alignment results.

Consider, for example, the following:
`bowtie2 -x refgen -U reads_1.fastq -S eg1.sam`

#### Are we done yet?

Not quite.
In order to properly hand over sequence alignment and mapping files to most
downstream programs (including GroopM), we must first convert our SAM files
to BAM files, sort them, and then index them.

For this, a sequence like the following will suffice:
```
samtools view -bS eg1.sam > eg1.bam
samtools sort -f eg1.bam eg1.sorted.bam
samtools index eg1.sorted.bam
```

The first step translates our SAM file into a more efficient binary format
known as BAM.
The next sorts entries in our BAM file by their starting nucleotide in our
overall assembly.
The last step creates a suitably named index file for our sorted BAM file
that can be used by programs like GroopM to efficiently access parts of our
(large) BAM file.
