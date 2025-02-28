# IMPACT

This project is no longer in development, but instead being replaced by a more sophisticated algorithm (IMPAQT) that is currently in development. This repository remains so that those referencing my Undergraduate Thesis have the code available. 

## Introduction

IMPACT (Identifies Multiple Peaks and Counts Transcripts) is a gene 
expression quantification method for TAGseq experiments developed 
by Bradley Jenner for his Undergraduate Honors Thesis at UC Davis. 
It operates on assumptions made about the distribution of sequencing reads
along the 3' UTR of a gene to cluster reads and assign their combined 
read count to the most appropriate gene. This method is particularly
useful in non-model organisms where 3' UTRs for most genes are poorly
annotated, resulting in massive data loss. It also can generate a GTF file
defining the boundaries and expression levels for each identified cluster. 
It is a C++ tool that can take advantage of multiple threads and is more or
less the developers first attempt at making more performant software. For
a more detailed explanation of the algorithm and validation methods, 
please read the associated paper submitted for the thesis, Jenner_Undergradaute_Thesis.pdf.


This tool is still very much in development. It relies partially
on the bamtools and seqan C++ libraries, although mainly for parsing
and manipulating bam and gtf files. Additionally, the threadsafe queue
was possible thanks to EmbeddedArtistry. This reliance will be revisited in
future versions along with implementing a more sopisiticated clustering 
algorithm to improve read groupings and their gene assignment.

For questions or comments, please contact
Bradley Jenner at <bnjenner@ucdavis.edu>

## Installation

0. Make sure cmake and make are installed on your machine.

1. Clone this repository and change into it.
```
git clone https://github.com/bnjenner/impact.git
cd impact
```

2. Create a build directory and change into it.
```
mkdir build
cd build
```

3. Compile
```
cmake ../
make
```

4. Add path to bash profile
```
echo "export PATH=$PATH:path/to/build_directory" >> ~/.bash_profile
source ~/.bash_profile
```
5. Give it a go! 

## Usage
```
   impact [input.sorted.bam] [annotation.gtf|annotation.gff] [options]

DESCRIPTION:
    
   Identifies expressed transcripts using clusters of mapped reads from TAGseq experiments.
   Generates a counts file written to stdout and optionally a GTF file of identified read clusters.

PARAMETERS:

    -h, --help
          Display the help message.

    -t, --threads INTEGER
          Number of processes for multithreading. Default: 1.

    -l, --library-type STRING
          Library type. Paired end is not recommended. Only used to check proper pairing. One of
          single and paired. Default: single.

    -s, --strandedness STRING
          Strandedness of library. One of forward and reverse. Default: forward.

    -n, --nonunique-alignments
          Count primary and secondary read alignments.

    -q, --mapq-min INTEGER
          Minimum mapping quality score to consider for counts. Default: 1.

    -f, --feature-tag STRING
          Name of feature tag. Default: exon.

    -i, --feature-id STRING
          ID of feature (use for GFFs). Default: gene_id.

    -o, --output-gtf STRING
          Output read cluster GTF file and specify name.

    --version
          Display version information.

```
