Duplex-Sequencing
=================
README

These programs are meant to be run in order, and result in the transformation of two input .FASTQ files containing data from two reads of an Ilumina sequencing run (or another next-generation sequencer) into a paired-end BAM file containing duplex concensus sequences (DCS).  They will also generate a file containing each different tag present, and how many times it occured, as well as a file, called extraConsensus.bam containing single strand concensus sequences (SSCSs) that, for one reason or another, didn't have a mate (This is possible even if all reads origionaly had mates if the mate had too few or too many reads).  

New in this version:
Ability to run paired-end files
Choose which reads to make consensuses with  (dual_map v. mono_map)
Filter for 'good' reads
Filter for 'good' SSCS's (not too many N's)
Filter for most common cigar string
Filter out reads with soft-clipping
Customizable file names
Bash script writing program

Required packages/programs
BWA (written with V 0.6.2)
Samtools (written with V 0.1.17)
Python (written with V 2.7.3)
Pysam (written with V 0.7.5)
BioPython (written with V 1.62)

Instructions: 

First, run PE_BASH_MAKER.py, then make the bash script (.sh file) exicutable using the command below.  

chmod +x PE_DCS_CALC.*.*.sh

Run the bash script with:

bash [scriptname.sh] 3>&1 1>&2 2>&3 | tee -a log.txt

This should run the rest of the process through to an output paired-end BAM file, copying the contents of stderr to a log file for documentation and reporting purposes.  It is strongly sugested that the final sorted BAM file undergo post-processing with picard-tools-1.70/AddOrReplaceReadGroups.jar and GATK/GenomeAnalysisTK.jar, before generating statistics.  Do not run the bash script in a folder containing pre-existing SAM files, as it will delete them.  

Inputs:
	read-1-raw-data.fq
	read-2-raw-data.fq

Outputs:
	SSCS-output-file.bam
		PE SSCSs

	DCS-output-file.bam
		PE DCSs

	tag-counts-file.tagcounts
		Number of members in each group, regardless of weather or not a SSCS was created from that group.

	SSCS-output-file_UP.bam
		SSCSs that didn't have a mate
	
	SSCS-output-file_NM.bam
		Reads that were not even considered
	
	SSCS-output-file_LCC.bam
		Reads that had a less common cigar string than the one used for SSCS generation

	DCS-output-file_UP.bam
		DCSs for which the switchtag resulted in too many N's.  
	

PE_BASH_MAKER.py

This program makes a shell script so that the user will not need to enter the commands for all the programs himself.  When using it, navigate to the folder with your data, with all the programs in a different folder, and enter a relative path to that folder from the folder with your data.  This way, the program will be able to auto-build the path to the programs.  

usage: PE_BASH_MAKER.py [-h] [--ref REF] [--r1src R1SRC] [--r2src R2SRC]
                        [--min MINMEM] [--max MAXMEM] [--cut CUTOFF]
                        [--Ncut NCUT] [--rlength RLENGTH]
                        [--read_type READ_TYPE] [--spacers SPACERS]
                        [--slength SLENGTH] [--blength BLENGTH]
                        [--plengths PLENGTHS]

arguments:
  -h, --help            show this help message and exit
  --ref REF             .FASTA file containing the reference genome
  --r1src R1SRC         .fq file containing the raw read1 data
  --r2src R2SRC         .fq file containing the raw read2 data
  --min MINMEM          Minimum members for SSCS consensus [3]
  --max MAXMEM          Maximum members for SSCS consensus[1000]
  --cut CUTOFF          Mimimum percent matching for base choice in SSCS
                        consensus [0.8]
  --Ncut NCUT           Maxumum percent N's allowed [0.1]
  --rlength RLENGTH     Length of a single read [80]
  --read_type READ_TYPE
                        Type of read. Options: 
							dual_map: both reads map propperly. Doesn't consider read pairs where only one maps. 
							mono_map: considers any read pair where one read maps. 
							[dual_map]

ConsensusMaker2.2.py

Consensus Maker
Version 2.2
By Brendan Kohrn and Scott Kennedy(1)
(1) Department of Pathology, University of Washington School of Medicine, Seattle, WA 98195
August 29, 2013


Written for Python 2.7.3
Required modules: Pysam, Samtools

This program is intended to be run on a paired-end BAM file, sorted by read position, with duplex tags in the header and constant read length.  It will output a paired-end BAM file with single strand consensus sequences (SSCSs), and a .tagcounts file which contains the different tags (on both strands) and how many times they occur, even if they are not used in SSCS generation, in order by read.  In addition, it will output a BAM file of SSCSs which are unpaired, either because one of the pair didn't match the criteria for allignment, or because of some other reason, and a BAM file of all unconsidered sequences in the original file.  Quality scores on the output BAM files are meaningless.  The file produced by this program is meant to continue on through the duplex maker.  

The program starts at the position of the first good read, determined by the type of read specified on startup.  It then goes through the file until it finds a new position, saving all reads as it goes.  When it finds a new position, it sends the saved reads to the consensus maker, one tag at a time, untill it runs out of tags.  Consensus sequences are saved until their mates come up, at which point both are written to the output BAM file, first read first.  After emptying the reads from the first position, it continues on through the origional file until it finds another new position, sends those reads to the consensus maker, and so on until the end of the file.  At the end of the file, any remaining reads are sent through the consensus maker, and any unpaired consensuses are written to extraConsensus.bam.  

In the future, the program may be able to autodetect read length.  

usage: ConsensusMaker2.2.py [-h] [--infile INFILE] [--tagfile TAGFILE]
                            [--outfile OUTFILE] [--rep_filt REP_FILT]
                            [--minmem MINMEM] [--maxmem MAXMEM]
                            [--cutoff CUTOFF] [--Ncutoff NCUTOFF]
                            [--readlength READ_LENGTH] [--read_type READ_TYPE]

arguments:
  -h, --help            show this help message and exit
  --infile INFILE       input BAM file
  --tagfile TAGFILE     output tagcounts file
  --outfile OUTFILE     output BAM file
  --rep_filt REP_FILT   Remove tags with homomeric runs of nucleotides of
                        length x [9]
  --minmem MINMEM       Minimum number of reads allowed to comprise a
                        consensus. [0]
  --maxmem MAXMEM       Maximum number of reads allowed to comprise a
                        consensus. [100]
  --cutoff CUTOFF       Percentage of nucleotides at a given position in a
                        read that must be identical in order for a consensus
                        to be called at that position. [0]
  --Ncutoff NCUTOFF     Maximum percentage of Ns allowed in a consensus [1]
  --readlength READ_LENGTH
                        Length of the input read that is being used. [80]
  --read_type READ_TYPE
                        Type of read. 
                        Options: 
							dual_map: both reads map propperly.  Doesn't consider read pairs where only one read maps. 
							mono_map: considers any read pair where one read maps. 


DuplexMaker2.2.py

DCS Maker
Version 2.2
By Brendan Kohrn and Scott Kennedy(1)
(1) Department of Pathology, University of Washington School of Medicine, Seattle, WA 98195	
August 23, 2013	

Written for\ Python 2.7.3
Required modules: Pysam, Samtools, BioPython

This program is intended to be run on a paired-end BAM file, sorted by read position, which has already been through the consensus maker.  It alligns SSCS's to their switchtag, and outputs a paired-end BAM file containing Duplex Consensus Sequences (DCS's) and a BAM file containing unpaired duplex consensus sequences.  

usage: DuplexMaker2.2.py [-h] [--infile INFILE] [--outfile OUTFILE]
                         [--Ncutoff NCUTOFF] [--readlength READ_LENGTH]

arguments:
  -h, --help            show this help message and exit
  --infile INFILE       input BAM file
  --outfile OUTFILE     output BAM file
  --Ncutoff NCUTOFF     Maximum percentage of Ns allowed in a consensus [1]
  --readlength READ_LENGTH
                        Length of the input read that is being used.  [80]
