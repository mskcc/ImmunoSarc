#!/bin/bash

fastq=$1
#jobname=$2
outdir=$2
stardir=$3
starindexdir=$4


SampleName=`echo $fastq | awk -F "/" '{print $NF}' | awk -F "_IGO" '{print $1}'`
echo -e $SampleName

bsub -cwd $outdir -J $SampleName -oo %J.o -eo %J.e -R "rusage[mem=40]" \
	${stardir}/source/STAR --genomeDir $starindexdir \
		--readFilesCommand zcat --readFilesIn $fastq/*fastq.gz --outFileNamePrefix $SampleName \
		--outSAMtype BAM SortedByCoordinate --runThreadN 20 \
		--outSAMattrRGline ID:$SampleName PL:Illumina PU:XXX LB:XXX SM:$SampleName CN:YYY \
		--outSAMunmapped Within \
    		--outFilterMultimapNmax 1 --outFilterMismatchNmax 3 \
    		--chimSegmentMin 10 --chimOutType WithinBAM SoftClip --chimJunctionOverhangMin 10 --chimScoreMin 1 --chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 \
		samtools index "$outdir/${SampleName}Aligned.sortedByCoord.out.bam"




