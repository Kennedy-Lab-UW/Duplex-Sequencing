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
samtools view -F 4 -b $1 > ${1/.aln.sort.bam/.filt.bam}

#----------------clipping final file---------------
echo 'Clipping Final File:'
java -jar $picard/AddOrReplaceReadGroups.jar INPUT=${1/.aln.sort.bam/.filt.bam} OUTPUT=${1/.aln.sort.bam/.filt.readgroups.bam} RGLB=UW RGPL=Illumina RGPU=ATATAT RGSM=default

samtools index ${1/.aln.sort.bam/.filt.readgroups.bam}

java -Xmx2g -jar $gaTK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $refGenome -I ${1/.aln.sort.bam/.filt.readgroups.bam} -o ${1/.aln.sort.bam/.filt.readgroups.intervals}

java -Xmx2g -jar $gaTK/GenomeAnalysisTK.jar -T IndelRealigner -R $refGenome -I ${1/.aln.sort.bam/.filt.readgroups.bam} -targetIntervals ${1/.aln.sort.bam/.filt.readgroups.intervals}  -o ${1/.aln.sort.bam/.filt.readgroups.realign.bam}

java -Xmx2g -jar $gaTK/GenomeAnalysisTK.jar -T ClipReads -I ${1/.aln.sort.bam/.filt.readgroups.realign.bam} -o ${1/.aln.sort.bam/.filt.readgroups.clipped.bam} -R $refGenome --cyclesToTrim "1-4,81-84" --clipRepresentation SOFTCLIP_BASES

#----------------generating stats from final file---------------
samtools mpileup -B -A -d 500000 -f $refGenome ${1/.aln.sort.bam/.filt.readgroups.clipped.bam} | tee ${1/.aln.sort.bam/.pileup} | python $progPath/CountMuts.py -d $3 -c $4 -C $5 -u | tee ${1%%.aln.sort.bam}.d${3}-c${4}-${5}.unique.countmuts

python $progPath/mut-position.py -i ${1/.aln.sort.bam/.pileup} -o ${1%%.aln.sort.bam}.d${3}-c${4}-${5}.mutpos -d $3 -c $4 -C $5
