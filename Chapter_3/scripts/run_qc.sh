#!/bin/bash

# run_qc.sh
# usage: bash run_qc.sh

FILES=*.fastq.gz
for file in $FILES
do
	filename=$(basename "$file")
	filename="${filename%.*}"
	echo "File on the loop: 	$filename"

# exclude samples that have already been processed by fastqc
#export PATH="/usr/local/opt/python/libexec/bin:$PATH"
#for f in /Volumes/Data_disk/bovis_lung_tnseq/*R1_001.fastq.gz
#do
#  file=`basename ${f%%.*}`
#  if [ ! -f ~/git/tn_seq/lung_tnseq/fastqc/${file}_fastqc.html ]
#  then
#     echo "processing sample ${file}"
     #call fastQC quality analysis
     fastqc ${file} -o ~/git/tn_seq/lung_tnseq/fastqc/gene_wiz2
#  fi
  echo -e "########################\n\n"
done


# Run MultiQC
# -f overwrites existing files, . runs with files in current directory, -o output directory
echo "Running MultiQC..."
export PATH="/usr/local/opt/python/libexec/bin:$PATH"
multiqc -f ~/git/tn_seq/lung_tnseq/fastqc/gene_wiz2
