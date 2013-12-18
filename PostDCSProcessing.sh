#!/bin/bash

# usage: bash /PATH/PostDCSProcessing.sh in_file ref_genome minDepth minClonality maxClonality

clear
#Path to PicardTools:
picard=~/Desktop/bioinformatics/programs/picard-tools-1.70

#Path to GATK:
gaTK=~/Desktop/bioinformatics/programs/GATK

#Path to reference genome:
refGenome=$2 #enter full path

progPath=${0%%PostDCSProcessing.sh}
echo 'Generating Analysis'
echo ' '
echo 'File is: ' $1
echo 'Reference Gemone is: '$refGenome
echo 'Minimum Depth: '$3
echo 'Minimum Clonality: '$4
echo 'Maximum Clonality: '$5
echo ' '

#----------------filter for maping reads---------------
echo 'Filtering:'
samtools view -F 4 -b $1 > ${1/.bam/.filt.bam}

#----------------clipping final file---------------
echo 'Clipping Final File:'
java -jar $picard/AddOrReplaceReadGroups.jar INPUT=${1/.bam/.filt.bam} OUTPUT=${1/.bam/.filt.readgroups.bam} RGLB=UW RGPL=Illumina RGPU=ATATAT RGSM=default

samtools index ${1/.bam/.filt.readgroups.bam}

java -Xmx2g -jar $gaTK/GenomeAnalysisTK.jar -T ClipReads -I ${1/.bam/.filt.readgroups.bam} -o ${1/.bam/.clipped.bam} -R $refGenome --cyclesToTrim "1-5,80-84" --clipRepresentation SOFTCLIP_BASES --fix_misencoded_quality_scores

#----------------generating stats from final file---------------
samtools mpileup -B -A -d 500000 -f $refGenome ${1/.bam/.clipped.bam} | tee ${1/.bam/.pileup} | python $progPath/CountMuts.py -d $3 -c $4 -C $5 -u | tee ${1%%.bam}.d${3}-c${4}-${5}.unique.countmuts

python $progPath/mut-position.py -i ${1/.bam/.pileup} -o ${1%%.bam}.d${3}-c${4}-${5}.mutpos -d $3 -c $4 -C $5
