#!/bin/bash -e
#SBATCH -c 4 --mem=10Gb --time 24:00:00 -J 12566_PILON_EV

module load Pilon/1.24-Java-15.0.2
module load BWA/0.7.18-GCC-12.3.0
module load SAMtools/1.19-GCC-12.3.0


mkdir -p Polished_Assemblies/

GENOMES=("23EVAB351" "23EV558" "23EV565" "23EV573" "23EV574" "23EV585" "23EV599" "23EV612" "23EV634" "23EV652" "23EV674")


for GENOME_PREFIX in "${GENOMES[@]}"; do
  echo "Processing genome $GENOME_PREFIX..."

GENOME="${GENOME_PREFIX}_contigs.fasta"
READS_1="Illumina_Assemblies/${GENOME_PREFIX}_R1.fq.gz" 
READS_2="Illumina_Assemblies/${GENOME_PREFIX}_R2.fq.gz"

bwa index $GENOME #index reference genome

bwa mem -t 8 $GENOME $READS_1 $READS_2 | samtools view -Sb - | samtools sort -o ${GENOME_PREFIX}.sorted.bam ##Align Illumina reads to genome 

samtools index ${GENOME_PREFIX}.sorted.bam #Index bam file

java -Xmx10G -jar $EBROOTPILON/pilon.jar \
 --genome $GENOME \
--fix all \
--changes \
--frags ${GENOME_PREFIX}.sorted.bam \
--threads 10 \
 --output Polished_Assemblies/${GENOME_PREFIX}_polished | tee Polished_Assemblies/round1_${GENOME_PREFIX}.pilon 

done 
