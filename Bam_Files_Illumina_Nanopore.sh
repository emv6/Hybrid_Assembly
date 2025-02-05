#!/bin/bash -e
#SBATCH -c 4 --mem=10Gb --time 24:00:00 -J 12566_PILON_EV

module load BWA/0.7.18-GCC-12.3.0
module load SAMtools/1.19-GCC-12.3.0
module load minimap2/2.28-GCC-12.3.0

GENOMES=("23EVAB351" "23EV558" "23EV565" "23EV573" "23EV574" "23EV585" "23EV599" "23EV612" "23EV634" "23EV652" "23EV674")


for GENOME_PREFIX in "${GENOMES[@]}"; do
  echo "Processing genome $GENOME_PREFIX..."

ILLUMINA_GENOME="${GENOME_PREFIX}_contigs.fasta"
READS_1="Illumina_Assemblies/${GENOME_PREFIX}_R1.fq.gz" 
READS_2="Illumina_Assemblies/${GENOME_PREFIX}_R2.fq.gz"
NANOPORE_READS="${GENOME_PREFIX}_Chopper.fastq.gz"
NANOPORE_GENOME="${GENOME_PREFIX}_NanoporeAssembly.fasta"

#Illumina Assembly
bwa index $ILLUMINA_GENOME #index reference genome
bwa mem -t 8 $ILLUMINA_GENOME $READS_1 $READS_2 | samtools view -Sb - | samtools sort -o ${GENOME_PREFIX}_Illumina_Illumina.sorted.bam ##Align Illumina reads to genome 
samtools depth -a ${GENOME_PREFIX}_Illumina_Illumina.sorted.bam | awk '{sum+=$3} END {print "Coverage = ",sum/NR}' #Illumina Coverage Depth

minimap2 -t 8 -ax map-ont $NANOPORE_GENOME $NANOPORE_READS | samtools view -Sb - | samtools sort -o ${GENOME_PREFIX}_nanopore.bam #Align nanopore reads to nanopore genome
samtools depth -a ${GENOME_PREFIX}_nanopore.bam | awk '{sum+=$3} END {print "Coverage = ",sum/NR}' #Nanopore Coverage Depth

done 
