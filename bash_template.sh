#!/bin/bash
# DS bash script
# Version 2.1
# 

# Step 1: Setup variables for run:
clear

# Set up error checking
# Stop on any error
set -e
# Stop on an error inside a pipeline
set -o pipefail
# Throw an error on calling an unassigned variable
set -u

#DEFAULTS
DSpath=''
alignRef=''
runIdentifier=''
read1in=seq1.fq
read2in=seq2.fq
iSize=-1
minMem=3
maxMem=1000
cutOff=0.7
nCutOff=1
readLength=0
barcodeLength=12
spacerLength=5
filtersSet='os'
readTypes='dpm'
repFilt=9
readOut=1000000
Ncores=18

#NONDEFAULTS

#FINAL_READ_LENGTH
readLength=$((readLength-barcodeLength-spacerLength))

#LOG_FILE_NAME
logFile=${runIdentifier}.log.txt

#Output folder
output_folder=${runIdentifier}

#Export all variables
export DSpath
export alignRef
export runIdentifier
export read1in
export read2in
export iSize
export minMem
export maxMem
export cutOff
export nCutOff
export readLength
export barcodeLength
export spacerLength
export filtersSet
export readTypes
export repFilt
export readOut
export Ncores
export output_folder

# Load required software into path using the Environment Modules Project (http://modules.sourceforge.net)
module load Python
module load BWA
module load SAMtools

# Print out options used to log file
touch $logFile
echo "Run identifier: " $runIdentifier | tee -a ${logFile}
echo "Program path: " $DSpath | tee -a ${logFile}
echo "Reference genome: " $alignRef | tee -a ${logFile}
echo "Barcode length: " $barcodeLength | tee -a ${logFile}
echo "Spacer length: " $spacerLength | tee -a ${logFile}
echo "Post-tag_to_header read length: " $readLength | tee -a ${logFile}
echo "Repetitive tag filter length: " $repFilt | tee -a ${logFile}
echo "Minimum family size: " $minMem | tee -a ${logFile}
echo "Maximum family size: " $maxMem | tee -a ${logFile}
echo "Consensus cutoff: " $cutOff | tee -a ${logFile}
echo "Consensus N cutoff: " $nCutOff | tee -a ${logFile}
echo "Read types: " $readTypes | tee -a ${logFile}
echo "Filters: " $filtersSet | tee -a ${logFile}
echo "" | tee -a ${logFile}

#  Step 2: Run tag_to_header.py on imput files

echo "Starting Run" | tee -a ${logFile}
echo "tag_to_header starting"  | tee -a ${logFile}
date | tee -a ${logFile}
echo "" | tee -a ${logFile}

python ${DSpath}/tag_to_header.py --infile1 $read1in --infile2 $read2in --outprefix ${runIdentifier} --tagstats --spacerlen ${spacerLength} --taglen ${barcodeLength}

# Step 3: Align sequences

echo "Aligning with BWA" | tee -a ${logFile}
date | tee -a ${logFile}

bwa aln -t ${Ncores} $alignRef ${runIdentifier}.seq1.smi.fq > ${runIdentifier}.seq1.aln
bwa aln -t ${Ncores} $alignRef ${runIdentifier}.seq2.smi.fq > ${runIdentifier}.seq2.aln
bwa sampe -s $alignRef ${runIdentifier}.seq1.aln ${runIdentifier}.seq2.aln ${runIdentifier}.seq1.smi.fq ${runIdentifier}.seq2.smi.fq > ${runIdentifier}.pe.sam

# Step 4: Sort aligned sequences
echo "Sorting aligned sequences" | tee -a ${logFile}
date | tee -a ${logFile}

samtools view -Sbu ${runIdentifier}.pe.sam | samtools sort - ${runIdentifier}.pe.sort

# Step 5: Run Consensus Maker
echo "Starting Consensus Maker" | tee -a ${logFile}
date | tee -a ${logFile}

python ${DSpath}/ConsensusMaker.py --infile ${runIdentifier}.pe.sort.bam --tagfile ${runIdentifier}.pe.tagcounts --outfile ${runIdentifier}.sscs.bam --minmem $minMem --maxmem $maxMem --readlength $readLength --cutoff $cutOff --Ncutoff $nCutOff --read_type $readTypes --filt $filtersSet --isize $iSize

# Step 6: Sort SSCSs
echo "Sorting SSCSs" | tee -a ${logFile}
date | tee -a ${logFile}

samtools view -bu ${runIdentifier}.sscs.bam | samtools sort - ${runIdentifier}.sscs.sort

# Step 7: Run Duplex Maker
echo "Starting Duplex Maker" | tee -a ${logFile}
date  | tee -a ${logFile}

python ${DSpath}/DuplexMaker.py --infile ${runIdentifier}.sscs.sort.bam --outfile ${runIdentifier}.dcs --Ncutoff $nCutOff --readlength $readLength

# Step 8: Align DCSs
echo "Aligning DCSs" | tee -a ${logFile}
date | tee -a ${logFile}

bwa aln -t ${Ncores}  $alignRef ${runIdentifier}.dcs.r1.fq > ${runIdentifier}.dcs.r1.aln
bwa aln -t ${Ncores} $alignRef ${runIdentifier}.dcs.r2.fq > ${runIdentifier}.dcs.r2.aln
bwa sampe -s $alignRef ${runIdentifier}.dcs.r1.aln ${runIdentifier}.dcs.r2.aln ${runIdentifier}.dcs.r1.fq ${runIdentifier}.dcs.r2.fq > ${runIdentifier}.dcs.sam

# Step 9: Sort aligned DCSs
echo "Sorting aligned DCSs" | tee -a ${logFile}
date | tee -a ${logFile}

samtools view -Sbu ${runIdentifier}.dcs.sam | samtools sort - ${runIdentifier}.dcs.aln.sort

# Step 10: Index sorted DCSs
echo "Indexing sorted DCSs" | tee -a ${logFile}
date | tee -a ${logFile}

samtools index ${runIdentifier}.dcs.aln.sort.bam

# Step 11: Clean up
echo "Finishing with run.. " $runIdentifier | tee -a ${logFile}
echo "Cleaning.." | tee -a ${logFile}
date | tee -a ${logFile} 
python ${DSpath}/clean.py --scripts_folder $(pwd) --output_folder ${output_folder} 
