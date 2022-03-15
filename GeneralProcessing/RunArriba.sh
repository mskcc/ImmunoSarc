#!/bin/bash

bampath=$1
outdir=$2
arribadir=$3
genomeFasta=$4
genomeGtf=$5


mkdir $outdir

bsub -cwd $outdir -oo %J.o -eo %J.e -R "rusage[mem=30]" -We 00:59 \
	"${arribadir}/arriba" \
    		-x $bampath \
    		-o fusions.tsv -O fusions.discarded.tsv \
    		-a $genomeFasta -g $genomeGtf \
		-b "${arribadir}/database/blacklist_hg19_hs37d5_GRCh37*" \
		-T -P

