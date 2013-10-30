#!/bin/bash
#Arguments
#in_file ref_genome minDepth minClonality maxClonality

clear
#Path to PicardTools:
picard=/Users/LoeblabM11/Desktop/bioinformatics/programs/picard-tools-1.70

#Path to GATK:
gaTK=/Users/LoeblabM11/Desktop/bioinformatics/programs/GATK

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
samtools mpileup -B -d 500000 -f $refGenome ${1/.aln.sort.bam/.filt.bam} | tee ${1/.aln.sort.bam/.filt.pileup} | python $progPath/CountMuts.py -d $3 -c $4 -C $5 > ${1%%.aln.sort.bam}.filt.pileup.d${3}-c${4}-${5}.countmuts
cat ${1/.aln.sort.bam/.filt.pileup} | python $progPath/CountMuts.py -d $3 -c $4 -C $5 -u > ${1%%.aln.sort.bam}.filt.pileup.d${3}-c${4}-${5}.unique.countmuts

samtools mpileup -B -d 500000 -f $refGenome ${1/.aln.sort.bam/.filt.readgroups.clipped.bam} | tee ${1/.aln.sort.bam/.filt.readgroups.clipped.bam.pileup} | python $progPath/CountMuts.py -d $3 -c $4 -C $5 > ${1%%.aln.sort.bam}.filt.readgroups.clipped.bam.pileup.d${3}-c${4}-${5}.countmuts
cat ${1/.aln.sort.bam/.filt.readgroups.clipped.bam.pileup} | python $progPath/CountMuts.py -d $3 -c $4 -C $5 -u > ${1%%.aln.sort.bam}.filt.readgroups.clipped.bam.pileup.d${3}-c${4}-${5}.unique.countmuts

cat ${1/.aln.sort.bam/.filt.readgroups.clipped.bam.pileup} | python $progPath/mut-position.py -d $3 -c $4 -C $5 > ${1%%.aln.sort.bam}.filt.readgroups.clipped.bam.pileup.d${3}-c${4}-${5}.mutpos
