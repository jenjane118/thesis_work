#!/bin/bash

#run_tpp.sh

#usage:
#conda activate tn_seq; bash ~/git/tn_seq/scripts/run_tpp.sh list_of_samples.txt

#read list of fastq files and analyse with tpp

FILE_LIST=$1
REF_FILE=~/tn_seq/data/Mbovis_AF2122_97.fasta
BWA_EXEC=~/anaconda3/envs/tnseq/bin/bwa
OUTPUT="lung_tpp"

for f in `cat $FILE_LIST`;
do
   echo "File in process: ${f}"
   sample=${f##*/}
   sample=${sample/_R1_001.fastq/}
   echo "Sample= $sample"
   #tpp -bwa $BWA_EXEC -ref $REF_FILE -reads1 ${f} -output "$OUTPUT"_"$sample"
   echo	"File complete= ${OUTPUT}_${sample}"
   echo "$OUTPUT"_"$sample"
#export PATH="/usr/local/opt/python/libexec/bin:$PATH"
done
