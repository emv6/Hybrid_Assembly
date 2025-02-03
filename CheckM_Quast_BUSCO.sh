#!/bin/bash -e
#SBATCH --cpus-per-task=10 --mem 10Gb  --time 1:00:00 -J 12566_EV_QUALITY

module purge
module load BUSCO/5.8.2-gimkl-2022a
module load CheckM/1.2.1-gimkl-2022a-Python-3.10.5
module load QUAST/5.2.0-gimkl-2022a

GENOMES=("23EVAB351" "23EVEV558" "23EVEV565" "23EVEV573" "23EVEV574" "23EVEV585" "23EVEV599" "23EV612" "23EV634" "23EV652" "23EV674")
OUTPUT=CHECKMOUTPUT/

checkm taxon_set genus "Staphylococcus" staph.ms 

for GENOME_PREFIX in "${GENOMES[@]}"; do
  echo "Processing genome $GENOME_PREFIX..."

busco -m genome -l staphylococcus_odb12 -i {GENOME_PREFIX}_contigs.fasta
checkm analyze -t 12 -x fasta staph.ms {GENOME_PREFIX}_contigs.fasta $OUTPUT
checkm qa -t 12 staph.ms $OUTPUT
quast.py {GENOME_PREFIX}_contigs.fasta

done
