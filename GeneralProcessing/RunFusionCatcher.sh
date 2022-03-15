#!/bin/bash

# $1 - path to fastq
# $2 - output directory
# $3 - job name

inputfastq=$1
outputdir=$2
jobname=$3
fc=$4

mkdir $2

bsub -cwd $2 -o log/ -e $jobname.error/ -R "rusage[mem=4]" -n 16 -W 48:00 -J $jobname \
        ${fc}/bin/fusioncatcher -p 16\
        -d ${fc}/data/current \
        -i $inputfastq \
        -o $outputdir \
