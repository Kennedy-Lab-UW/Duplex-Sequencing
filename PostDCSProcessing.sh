#!/bin/bash

# usage: bash /PATH/PostDCSProcessing.sh in_file ref_genome minDepth minClonality maxClonality

clear
#Path to PicardTools:
picard=/PATH/picard-tools-1.70

#Path to GATK:
gaTK=/PATH/GATK

#Path to reference genome:
refGenome=$2 #enter full path

progPath=${0%%PostDCSProcessing.sh}
echo 'Generating Analysis'
echo
echo 'File is: ' $1
echo 'Reference Gemone is: '$refGenome
echo 'Minimum Depth: '$3
echo 'Minimum Clonality: '$4
echo 'Maximum Clonality: '$5
echo

#----------------filter for maping reads---------------
echo 'Filtering:'
samtools view -F 4 -b $1 > ${1/.aln.sort.bam/.filt.bam}

#----------------clipping final file---------------
echo 'Clipping Final File:'
java -jar $picard/AddOrReplaceReadGroups.jar INPUT=${1/.aln.sort.bam/.filt.bam} OUTPUT=${1/.aln.sort.bam/.filt.readgroups.bam} RGLB=UW RGPL=Illumina RGPU=ATATAT RGSM=default

java -Xmx2g -jar $gaTK/GenomeAnalysisTK.jar -T ClipReads -I ${1/.aln.sort.bam/.filt.readgroups.bam} -o ${1/.aln.sort.bam/.filt.readgroups.clipped.bam} -R $refGenome --cyclesToTrim "1-5,71-80" --clipRepresentation HARDCLIP_BASES --fix_misencoded_quality_scores --unsafe

#----------------generating stats from final file---------------
samtools mpileup -B -d 500000 -f $refGenome ${1/.aln.sort.bam/.filt.bam} | tee ${1/.aln.sort.bam/.filt.pileup} | python $progPath/CountMuts.py -o ${1%%.aln.sort.bam}.filt.pileup.d${3}-c${4}-${5}.countmuts -d $3 -c $4 -C $5
python $progPath/CountMuts.py -i ${1/.aln.sort.bam/.filt.pileup} -o ${1%%.aln.sort.bam}.filt.pileup.d${3}-c${4}-${5}.unique.countmuts -d $3 -c $4 -C $5 -u

samtools mpileup -B -d 500000 -f $refGenome ${1/.aln.sort.bam/.filt.readgroups.clipped.bam} | tee ${1/.aln.sort.bam/.filt.readgroups.clipped.bam.pileup} | python $progPath/CountMuts.py -o ${1%%.aln.sort.bam}.filt.readgroups.clipped.bam.pileup.d${3}-c${4}-${5}.countmuts -d $3 -c $4 -C $5
python $progPath/CountMuts.py -i ${1/.aln.sort.bam/.filt.readgroups.clipped.bam.pileup} -o ${1%%.aln.sort.bam}.filt.readgroups.clipped.bam.pileup.d${3}-c${4}-${5}.unique.countmuts -d $3 -c $4 -C $5 -u

python $progPath/mut-position.py -i ${1/.aln.sort.bam/.filt.readgroups.clipped.bam.pileup} -o ${1%%.aln.sort.bam}.filt.readgroups.clipped.bam.pileup.d${3}-c${4}-${5}.mutpos -d $3 -c $4 -C $5