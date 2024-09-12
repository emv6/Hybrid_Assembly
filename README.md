# HybridAssembly

![bash](https://img.shields.io/badge/language-bash-green)

All bioinformatic analysis was conducted on the New Zealand eScience Infrastructure [NeSI](https://github.com/nesi). FastQC used for QC of the Illumina sequencing of the 23 *S. aureus* samples.

## FastQC 
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) was used to check for QC of the samples for adaptor content and sequence quality

```#!/bin/bash -e
#SBATCH --cpus-per-task=8 --mem 50Gb --time 1:00:00 -J FASTQC_EV

module load FastQC/0.11.9

FASTQC = /pathtorawreads
OUTPUT_DIR = /processedreadsdirectory/FASTQC/

mkdir -p $OUTPUT_DIR

for file in $FASTQC/*.fq.gz;
do
if [-f $FILE ];
then echo "Running FastQC on $FILE"
fastqc $FILE -o $OUTPUTDIR
else
echo "No FASTQ files found in $FASTQ"
fi
done
echo "FastQC analysis completed for all samples"
```

## Assembly method of 23 *S. aureus* isolates 

#### Guppy
[Guppy](https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol/v/gpb_2003_v1_revax_14dec2018/guppy-software-overview) was chosen as the basecaller specifying the super accurate model. 
#### Filtlong
[Filtlong](https://github.com/rrwick/Filtlong) was used to remove reads shorter than 1kbp and exclude the worst 5% of the reads. 
#### Porechop
[Porechop](https://github.com/rrwick/Porechop) was used to remove any adaptors 
#### NanoStat
[NanoStat](https://github.com/wdecoster/nanostat) was used to perform quality checks to determine the mean read quality of each genome
#### Chopper
[Chopper](https://github.com/wdecoster/chopper) was used to ensure only reads above the mean read quality are kept. 
#### Flye 
[Flye](https://github.com/mikolmogorov/Flye) was used as the *de novo* assembler to assemble each genome
#### Medaka 
[Medaka](https://github.com/nanoporetech/medaka) was used as the polishing tool for each genome
#### BWA and SAMtools 
[BWA](https://github.com/lh3/bwa) and [SAMtools](https://samtools.sourceforge.net/) were used to iteratively map the Illumina reads of each genome to the nanopore assembly. 
#### Pilon
[Pilon](https://github.com/broadinstitute/pilon) is run to polish the final Hybrid Assembly. Three rounds of Pilon are run. Therefore Pilon is run three times and the bwamem_pilon is run twice. The final fasta file called pilon_round3.fasta is the polished hybrid genome and is renamed accordingly. 
#### CheckM and Quast 
[CheckM](https://github.com/Ecogenomics/CheckM) and [Quast](https://github.com/ablab/quast) was run on each genome to check the assembly is complete 
