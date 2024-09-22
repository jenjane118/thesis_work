#!/bin/bash
input=$1
while IFS= read -r line

do
	filename=$(basename "$line")
	extension="${filename##*.}"
	filename="${filename%.*}"

	echo "file on the loop:		$line"

	transit resampling ~/git/tn_seq/lung_tnseq/final_sample_wigs/perm_lung_tpp_MbA27.wig ${line} ~/git/tn_seq/ncrna_mbovis.prot_table resamp_ncrna_${filename}.txt

	echo "output:                    resamp_ncrna_${filename}.txt"

done <"$input"
