Guppy
```
#!/bin/bash -e

#SBATCH -J guppy_gpu --gpus-per-node P100:1 --mem 6G --cpus-per-task 8 --time 48:00:00 --output slurmout_%A_%a.out --error slurmout_%A_%a.err --array=0-20

module purge
module load ont-guppy-gpu/6.4.6

BASE_DIR=/Guppy/
DIRS=(${BASE_DIR}/K*/)

echo $DIRS

INPUT=${DIRS[$SLURM_ARRAY_TASK_ID]}
OUTPUT="${INPUT}GuppySuperAccurate"

echo $OUTPUT
mkdir -p $OUTPUT

guppy_basecaller -i "$INPUT" -s "$OUTPUT" -c /opt/nesi/CS400_centos7_bdw/ont-guppy-gpu/6.4.6/data/dna_r9.4.1_450bps_sup.cfg --device auto --recursive --detect_mid_strand_adapter

#Genome is complete
echo "Processing of $INPUT complete"```
```
cat $OUTPUT/*.fastq > $filename.fastq #change filename to isolate name
```

Filtlong 
```
module load Filtlong/0.2.0
filtlong --min_length 1000 --keep_percent 95 $filename_Guppy.fastq| gzip > $filename_Guppy_Filtlong.fastq.gz
```

Porechop
```
#!/bin/bash -e
#SBATCH -c 15 --mem=16Gb --time 00:120:00 -J Porechop_EV

module purge
module load Porechop/0.2.4-gimkl-2020a-Python-3.8.2

for input in *_Guppy_Filtlong.fastq.gz
do
base=$(basename ${input} _Guppy_Filtlong.fastq.gz)
echo "working with file $input"
echo "basename is $base"

OUTPUT=StaphA_${base}_Guppy_Filtlong_Porechop.fastq.gz
porechop -i $input -o $OUTPUT --discard_middle

done
```
NanoStat 
```
module load NanoStat/1.5.0-gimkl-2020a-Python-3.8.2
NanoStat --fastq StaphA_$filename_Guppy_Filtlong_Porechop.fastq.gz
```
Chopper
```
module load chopper/0.2.0-GCC-11.3.0
gunzip -c StaphA_$filename_Guppy_Filtlong_Porechop.fastq.gz | chopper -q $quality | gzip > StaphA_$filename_Guppy_Filtlong_Porechop_Chopper.fastq.gz
#The q is the mean read quality as specified from NanoStat - $quality is altered to the output from the mean read quality from NanoStat.
```
Flye
```
#!/bin/bash -e
#SBATCH --nodes 1 --cpus-per-task 1 --ntasks 10 --mem=50G -J Flye_EV --time=72:00:00 --hint=nomultithread

module load Flye/2.9.1-gimkl-2022a-Python-3.10.5
INPUT=StaphA_$filename_Guppy_Filtlong_Porechop_Chopper.fastq.gz
OUTPUT=$filename_Flye/
flye --nano-hq $INPUT --genome-size 2.8m -o $OUTPUT -t 10 -i 3 --asm-coverage 50
```
Medaka
```
#!/bin/bash -e
#SBATCH --nodes 1 --cpus-per-task 1 --job-name medaka_Flye --mem=20G --time=24:00:00 --output=%x_%j.out --error=%x_%j.err --hint=nomultithread

module load medaka/1.6.0-Miniconda3-4.12.0

READS=StaphA_$filename_Guppy_Filtlong_Porechop_Chopper.fastq.gz
CONTIG_FILE=$filename_Flye/assembly.fasta

medaka_consensus -i $READS -d $CONTIG_FILE -o medaka_$filename/ -t 2 -m r941_min_sup_g507
```
```
#Rename medaka polished output
mv medakaoutput.fasta $filename_Medaka_Polish.fasta
```
BWA & SAMtools
```
#!/bin/bash -e
#SBATCH -c 10 --mem=8Gb --time 00:10:00 -J bwamem_EV

module load BWA/0.7.17-GCC-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0

REF=23EV612_Medaka_Polish.fasta
bwa mem -t 6 -R"@RG\tID:23EV612\tPL:ILLUMINA\tSM:LIC_$filename_2019" $REF $filename_R1.fq.gz $filename_R2.fq.gz | samtools view - -Sb | samtools sort - -@10 -o $filename_Medaka_Illumina.bam
```
```
#The output bam is indexed using SAMtools 
module load SAMtools/1.16.1-GCC-11.3.0
samtools index $filename_Medaka_Illumina.bam
```
Pilon & BWAmem_Pilon
```
#!/bin/bash -e
#SBATCH -c 4 --mem=10Gb --time 00:10:00 -J PILON_EV

mkdir -p $filename_Pilon/
OUTPUT=$filename_Pilon/ 

module load Pilon/1.24-Java-15.0.2

GENOME=$filename_Medaka_Polish.fasta
BAM=$filename_Medaka_Illumina.bam 

java -Xmx10G -jar $EBROOTPILON/pilon.jar --genome $GENOME --fix all --changes --frags $BAM --threads 10 --output $OUTPUT/pilon_round1 | tee $OUTPUT/round1.pilon
```
```
#The output fasta file is indexed using BWA
module load BWA/0.7.17-GCC-11.3.0
samtools index pilon_round1.fasta 
```
```
#!/bin/bash -e
#SBATCH -c 5 --mem=8Gb --time 00:10:00 -J bwamempilon_EV

module load BWA/0.7.17-GCC-9.2.0
module load SAMtools/1.16.1-GCC-11.3.0

INPUT=pilon_round1.fasta
ILLUMINAR1=$filename_R1.fq.gz
ILLUMINAR2=$filename_R2.fq.gz
bwa mem -t 5 $INPUT $ILLUMINAR1 $ILLUMINAR2 | samtools view - -Sb | samtools sort - -@14 -o Pilon1_$filename.bam
```bash
The output bam is indexed using SAMtools 
module load SAMtools/1.16.1-GCC-11.3.0
samtools index Pilon1_$filename.bam
```
#Pilon was run three times and bwa_mem_pilon was run twice. 

