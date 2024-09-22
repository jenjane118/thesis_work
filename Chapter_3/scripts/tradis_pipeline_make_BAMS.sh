#!/bin/bash
# from Dong Xia's tradis pipeline

# usage: bash tradis_pipeline_make_BAMS.sh 

ref=/d/in16/u/sj003/refseqs/Mbovis_AF2122_97.fasta
tag=GTCTAGAGACCGGGGACTTATCAGCCAACCTGTTA


# indexes database sequences in FASTA format
bwa index ${ref}

while read i;
do
#gunzip ${i}.fastq.gz
# creates fastq file containing reads that match supplied tag
filter_tradis_tags -f ${i}.fastq -t ${tag} -o ${i}.removed0mis.fastq
# creates fastq file containing reads with supplied tag removed from seqs
remove_tradis_tags -f ${i}.removed0mis.fastq -t ${tag} -o ${i}.removed.fastq

# aln is old way: finds SA coordinates of input reads (what dong uses)
 bwa aln ${ref} ${i}.removed.fastq > ${i}.removed.fastq.sai
# aligns in SAM format
bwa samse ${ref} ${i}.removed.fastq.sai ${i}.removed.fastq > ${i}.sam

# bwa mem for 70bp-1Mbp query seqs with BWA-MEM algorithm using maximal exact matches and extending with SW algorithm (should have more reads) local alignment
# this is what TPP (DeJesus) uses
#bwa mem ${ref} ${i}.removed.fastq > ${i}.sam
# file format conversion
samtools view -bS ${i}.sam > ${i}.bam
# sort alignments by leftmost coordinates
samtools sort ${i}.bam > ${i}.sort.bam
# index sorted bam file for fast random access (output is aln.bam)
samtools index ${i}.sort.bam
# remove redundant files
rm ${i}.sam ${i}.bam ${i}.removed0mis.fastq ${i}.removed.fastq
# rename unmapped reads, sort, convert to fastq format
samtools view -b -f 4 ${i}.sort.bam > ${i}_unmapped.bam
samtools sort ${i}_unmapped.bam > ${i}_unmappedsort.bam
# use bedtools to convert bam files to fastq
bamToFastq -i ${i}_unmappedsort.bam -fq ${i}_unmapped.fastq
# rename mapped reads, sort, convert to fastq format
samtools view -b -F 4 ${i}.sort.bam > ${i}_mapped.bam
samtools sort ${i}_mapped.bam > ${i}_mappedsort.bam
bamToFastq -i ${i}_mappedsort.bam -fq ${i}_mapped.fastq
# remove redundant files
rm ${i}_unmapped.bam ${i}_unmappedsort.bam ${i}_mapped.bam ${i}_mappedsort.bam
done < names.txt

mkdir $(date +"%Y_%m_%d")_BAMs; mv *.bam *.bai $_
mkdir $(date +"%Y_%m_%d")_Mapped_read; mv *_mapped.fastq $_ 
mkdir $(date +"%Y_%m_%d")_Unmapped_reads; mv *_unmapped.fastq $_

