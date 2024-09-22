#!/bin/bash
# Runs fastp in PE mode for all .fastq files in a directory
# Run as:
# ./run_fastp.sh directory_of_samples 

timestamp=`date "+%Y%m%d-%H%M%S"`
logfile="run_$timestamp.log"
exec > $logfile 2>&1  #all output will be logged to logfile

FASTP_EXEC=~/anaconda3/envs/tnseq/bin/fastp
DIR=$1
ADAPTER_FILE=~/git/tn_seq/adapters_all.fa

echo "Running fastp using executable: $FASTP_EXEC"

for file in `ls $DIR/*_R1_*.fastq`; 
do
  sample=${file/$DIR\//}
  sample=${sample/_R1_001.fastq/}
  echo "Sample= $sample"
  $FASTP_EXEC -i  "$DIR/$sample"_R1_001.fastq  \
              -I  "$DIR/$sample"_R2_001.fastq \
         -o "$sample"_R1_001_trimmed.fastq \
         -O "$sample"_R2_001_trimmed.fastq \
         --adapter_fasta $ADAPTER_FILE \
         -l 30 --trim_poly_x \
         -h "$sample"_fastp.html -j "$sample"_fastp.json  
done


