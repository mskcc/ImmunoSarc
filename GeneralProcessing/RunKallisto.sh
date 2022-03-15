#!/bin/bash

# $1 - path to fastq R1 (.fq.qz)
# $2 - path to fastq R2 (.fq.qz)

Read1=$1
Read2=$2
ensembl_index=$3
kallisto=$4
outfold=$5
Sample_name=`echo $1 | awk -F "/" '{print $NF}' | awk -F "_IGO" '{print $1}'`
echo -e $Sample_name
out_dir=$5/$Sample_name


mkdir -p $5/$Sample_name
bsub -e $out_dir -n 5 -R rusage[mem=15] -We 1:59 $kallisto quant -i $ensembl_index -o $out_dir --bias -b 100 -t 5 --rf-stranded --fusion $Read1 $Read2

