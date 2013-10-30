Duplex-Sequencing
=================
README

Duplex Sequencing software package
Version 2.0
October 28, 2013
Programs by Scott Kennedy(1), Brendan Kohrn, and Mike Schmitt(1)
(1) Department of Pathology, University of Washington School of Medicine, Seattle, WA 98195
Several steps are based on prior work by Joe Hiatt

1. Glossery
    Duplex Consensus Sequence (DCS):
        A construct created by comparing two SSCSs.  Quality scores and cigar strings attached to DCS sequences are meaningless, though cigar strings regain meaning after reallignment.  
    Duplex tag:
        A random sequence of nucleotides that enables the identification of sequences resulting from the same starting molecule.  
    Family:  
        A group of reads that shares the same duplex tag. 
    Read: 
        A DNA sequence which has not been compressed by ConsensusMaker.py.  A raw read has not yet been modified by tag_to_header.py, while an SMI read has.
    Single Stranded Consensus Sequence (SSCS):
        A construct created by comparing multiple reads and deciding ambiguities by simple majority.  SSCSs are created by ConsensusMaker.py.  Quality scores attached to SSCSs are meaningless. although the cigar strings still have meaning.  

2. Summary of process
    These programs are meant to be run in order, and result in the transformation of two input .FASTQ files containing data from two reads of an Ilumina sequencing run (or another next-generation sequencer) into a paired-end BAM file containing DCSs.  They will also generate a file containing each different tag present, and how many times it occured, as well as a file, called extraConsensus.bam containing SSCSs that, for one reason or another, didn't have a mate (This is possible even if all reads origionaly had mates if the mate had too few or too many reads).  

3. Changes in this version from the last version (1.21)
    Ability to run paired-end files
    Bash script writing program
    Choose which reads to make consensuses with  (dual_map v. mono_map)
    Customizable file names
    Filter for 'good' reads
    Filter for 'good' SSCS's (not too many N's)
    Filter for most common cigar string
    Filter out reads with soft-clipping

4. Dependencies
    The following programs and packages must be installed on your computer.

    BWA (written with V 0.6.2)
    Samtools (written with V 0.1.17)
    Python (written with V 2.7.3)
    Pysam (written with V 0.7.5)
    BioPython (written with V 1.62)

5. Inputs:
	read-1-raw-data.fq
	read-2-raw-data.fq

6. Usage

    Create a folder with both your fastq files in it.  
    Run PE_BASH_MAKER.py, making sure to imput the correct read length (option --rlength), using the syntax shown below. Although it is recomended that all non-optional inputs be provided, the only inputs that are truely required are --ref, --r1src, and --r2src.  

           PE_BASH_MAKER.py [-h] [--ref REF] [--r1src R1SRC] [--r2src R2SRC]
                            [--min MINMEM] [--max MAXMEM] [--cut CUTOFF]
                            [--Ncut NCUT] [--rlength RLENGTH] [--blength BLENGTH]
                            [--slength SLENGTH] [--progInd PROGIND]
                            [--read_type READ_TYPE] [--isize ISIZE] [--absolute]
                            [--parallel]

    Arguments:
      -h, --help            show this help message and exit
      --ref REF             .FASTA file containing the reference genome
      --r1src R1SRC         .fq file containing the raw read1 data
      --r2src R2SRC         .fq file containing the raw read2 data
      --min MINMEM          Minimum members for SSCS consensus [3]
      --max MAXMEM          Maximum members for SSCS consensus [1000]
      --cut CUTOFF          Mimimum percent matching for base choice in SSCS
                            consensus [0.8]
      --Ncut NCUT           Maxumum percent N's allowed [0.1]
      --rlength RLENGTH     Length of a single read [85]
      --blength BLENGTH     length of the barcode sequence on a unprocessed single
                            read. [12]
      --slength SLENGTH     length of the spacer sequence in a unprocessed single
                            read.
      --progInd PROGIND     how often you want to be told what a program is doing
                            [1000000]
      --read_type READ_TYPE
                            Type of read. Options: dual_map: both reads map
                            properly. Doesn't consider read pairs where only one
                            read maps. mono_map: considers any read pair where one
                            read maps. [mono_map]
      --isize ISIZE         Optional: Maximum distance between read pairs [-1]
      --absolute            Optional: Treat the program path as an absolute path
      --parallel            Optional: Perform the alignments of both reads in
                            parallel

    Run the bash script from the command line with:

    bash scriptname.sh 3>&1 1>&2 2>&3 | tee -a log.txt

    where scriptname.sh is the newly-created bash script (.sh file), and log.txt is the desired name of your log file.  This should run the rest of the process through to an output paired-end BAM file, copying the contents of stderr to a log file for documentation and reporting purposes.  

    Note: Do not run the bash script in a folder containing pre-existing SAM files, as it will delete them.  

    It is strongly sugested that the final sorted BAM file undergo post-processing with picard-tools-1.70/AddOrReplaceReadGroups.jar and GATK/GenomeAnalysisTK.jar, before generating statistics.  

7. Data Outputs:
    * indicates a custome string representing one of the input files.  

	BAM file containing position-sorted paired-end reads: 
                                                PE.*.*.bam
    BAM file containing paired-end SSCSs:       SSCS.*.*.bam
    BAM file containing unpaired SSCSs:         SSCS.*.*_UP.bam
    BAM file containing non-mapping or otherwise bad reads: 
                                                SSCS.*.*_NM.bam
    BAM file containing good reads with less common cigar scores:
                                                SSCS.*.*_LCC.bam
    tagcounts file:                             PE.*.*.tagcounts
    BAM file containing paired-end DCSs:        DCS.*.*.bam
    BAM file containing unpaired DCSs:          DCS.*.*_UP.bam

    --note that these SSCS reads are all aligned relative to the reference genome, and have thus been reverse-complemented when necessary by the aligner. Thus these SSCS reads do NOT inform whether there is a strand bias due to DNA damage. Doing so requires looking at forward-mapping and reverse-mapping reads separately after the initial alignment. We intend to automate this type of analysis in a future version of our software.


8. Live Outputs

    The file Duplex-Process-Numbers.txt describes the number of reads in each file and the live outputs from each step.  

9. Analysis
    --While the main pipeline does no analysis, there are a number of options.  A shell script to perform analysis based on mutation frequency is provided, and approximates the analysis found in version 1.21.  If this script is to be used,  GATK and Picard Tools must be present on the computer, and the paths to both programs must be set propperly within the script.  Once this is done, the script can be run from any folder, using
        
        bash /PATH/PostDCSProcessing.sh DCS.*.*.aln.sort.bam /REFPATH/ref_genome.fasta minDepth minClonality maxClonality 2>&1 | tee -a log.txt
    
    where minDepth is an integer, minClonality is a decimal and maxClonality is a decimal.  

    --We have noticed that alignment errors at the ends of reads can result in false mutations. To eliminate these, we hard-clip the first and last 5 nt of each read after alignment: DCS.*.*.readgroups.clipped.bam

    --text file listing overall mutation frequencies: these are the files having extension .countmuts
    
    --text files listing unique mutation frequencies: these are the files having extension .unique.countmuts
    
    --text file listing position-specific mutation frequencies: files having extension .pileup.mutpos

    --the mutpos file is tab-delimited. Output is: 
        reference name, reference base, position number, depth, number of total mutations (excluding indels), number of mutations to T, C, G, A, insertions, deletions

10. Quality Control
    It is highly recommended to calculate read count statistics from each run for troubleshooting purposes. This can be done using a combination of samtools view and grep.
    
    You will want to consider the following numbers. If you lose a lot of data at a single step, you can then troubleshoot that step. For example if you lose a lot of reads going from 'mapped' to 'SSCS', the DNA is probably over-duplicated. Consider re-prepping the DNA using a larger amount of input into the PCR.

11. Debugging

    Should an error occur, it is most likely to occur durring the tag_to_header step (step 1), or durring analysis.  Problems in tag_to_header are likely to take the form of all reads comming out listed as bad reads.  This probably means that the file has gotten one line out of sync with the expected fastq format.  To locate the problem, you can re-run tag_to_header with a smaller --read_out parameter.  

    Should the program be stopped durring any step, it can be restarted by opening the bash script (.sh file) with your favorite text editor, commenting out the lines that have already run, and repeating the run command.  It is recomended that an entry be made in the log manually should something happen durring the run.  

    If an error occurs durring the analysis steps, consult the man page for the relevent programs.  

12. Details of the Individual Programs, in order of use (Advanced Users).  

    This information is also found at the top of each program, respectivly.  

    PE_BASH_MAKER.py
        PE Bash Maker V 1.0
        by Brendan Kohrn
        
        Write a bash script to run the process.  
        
        This program makes a shell script so that the user will not need to 
        enter the commands for all the programs himself.  When using it, 
        navigate to the folder with your data, with all the programs in a 
        different folder.  

        usage: PE_BASH_MAKER.py [-h] [--ref REF] [--r1src R1SRC] [--r2src R2SRC]
                                [--min MINMEM] [--max MAXMEM] [--cut CUTOFF]
                                [--Ncut NCUT] [--rlength RLENGTH] [--blength BLENGTH]
                                [--slength SLENGTH] [--progInd PROGIND]
                                [--read_type READ_TYPE] [--isize ISIZE] [--absolute]
                                [--parallel]

        Arguments:
          -h, --help            show this help message and exit
          --ref REF             .FASTA file containing the reference genome
          --r1src R1SRC         .fq file containing the raw read1 data
          --r2src R2SRC         .fq file containing the raw read2 data
          --min MINMEM          Minimum members for SSCS consensus [3]
          --max MAXMEM          Maximum members for SSCS consensus [1000]
          --cut CUTOFF          Mimimum percent matching for base choice in SSCS
                                consensus [0.8]
          --Ncut NCUT           Maxumum percent N's allowed [0.1]
          --rlength RLENGTH     Length of a single read [85]
          --blength BLENGTH     Length of the barcode sequence on a unprocessed single
                                read. [12]
          --slength SLENGTH     length of the spacer sequence in a unprocessed single
                                read.
          --progInd PROGIND     How often you want to be told what a program is doing
                                [1000000]
          --read_type READ_TYPE
                                Type of read. Options: dual_map: both reads map
                                properly. Doesn't consider read pairs where only one
                                read maps. mono_map: considers any read pair where one
                                read maps. [mono_map]
          --isize ISIZE         Optional: Maximum distance between read pairs [-1]
          --absolute            Optional: Treat the program path as an absolute path
          --parallel            Optional: Perform the alignments of both reads in
                                parallelOptional: Perform the alignments of both reads in parallel.  This is faster but requires more memory (minimum 16 GB recommended). 


    tag_to_header.py
        Tag To Header
        Version 2.0
        By Joe Hiatt, Scott Kennedy, Brendan Kohrn and Mike Schmitt
        October 23, 2013

        Isolate duplex tags, move them from within the sequenced read to the header region, and remove the spacer region.  

        usage: tag_to_header.py [-h] [--infile1 INFILE1] [--infile2 INFILE2]
                                [--outfile1 OUTFILE1] [--outfile2 OUTFILE2]
                                [--barcode_length BLENGTH] [--spacer_length SLENGTH]
                                [--read_out ROUT] [--adapter ADAPTERSEQ]

        optional arguments:
          -h, --help            show this help message and exit
          --infile1 INFILE1     First input raw fastq file.
          --infile2 INFILE2     Second input raw fastq file.
          --outfile1 OUTFILE1   Output file for first fastq reads.
          --outfile2 OUTFILE2   Output file for second fastq reads.
          --barcode_length BLENGTH
                                Length of the duplex tag sequence. [12]
          --spacer_length SLENGTH
                                Length of the spacer sequences used. [12]
          --read_out ROUT
          --adapter ADAPTERSEQ  Optional: Spacer sequence for filtering on the
                                presence of the spacer. This could be thrown off by
                                low quality scores.

    ConsensusMaker.py
        Consensus Maker
        Version 2.0
        By Brendan Kohrn and Scott Kennedy(1)
        (1) Department of Pathology, University of Washington School of Medicine, Seattle, WA 98195
        Based on work by Mike Schmitt and Joe Hiatt
        October 28, 2013

        Written for Python 2.7.3
        Required modules: Pysam, Samtools

        Inputs: 
            A position-sorted paired-end BAM file containing reads with a duplex tag in the header.  

        Outputs:
            1: A paired-end BAM file containing SSCSs
            2: A single-end BAM file containing unpaired SSCSs
            3: A single-end BAM file containing reads with less common cigar strings
            4: A single-end BAM file containing non-mapping reads
            5: A tagcounts file
            
            Note that quality scores in outputs 1, 2, and 3 are just space fillers and do not signify anything about the quality of the sequence.  

        The program starts at the position of the first good read, determined by the type of read specified on startup.  It then goes through the file until it finds a new position, saving all reads as it goes.  When it finds a new position, it sends the saved reads to the consensus maker, one tag at a time, untill it runs out of tags.  Consensus sequences are saved until their mates come up, at which point both are written to the output BAM file, read 1 first.  After making consensuses with the reads from the first position, it continues on through the origional file until it finds another new position, sends those reads to the consensus maker, and so on until the end of the file.  At the end of the file, any remaining reads are sent through the consensus maker, and any unpaired consensuses are written to a file ending in _UP.bam.  

        In the future, this program may be able to autodetect read length.  

        usage: ConsensusMaker.py [-h] [--infile INFILE] [--tagfile TAGFILE]
                                 [--outfile OUTFILE] [--rep_filt REP_FILT]
                                 [--minmem MINMEM] [--maxmem MAXMEM] [--cutoff CUTOFF]
                                 [--Ncutoff NCUTOFF] [--readlength READ_LENGTH]
                                 [--read_type READ_TYPE] [--isize ISIZE]
                                 [--read_out ROUT]

        optional arguments:
          -h, --help            show this help message and exit
          --infile INFILE       input BAM file
          --tagfile TAGFILE     output tagcounts file
          --outfile OUTFILE     output BAM file
          --rep_filt REP_FILT   Remove tags with homomeric runs of nucleotides of
                                length x. [9]
          --minmem MINMEM       Minimum number of reads allowed to comprise a
                                consensus. [3]
          --maxmem MAXMEM       Maximum number of reads allowed to comprise a
                                consensus. [1000]
          --cutoff CUTOFF       Percentage of nucleotides at a given position in a
                                read that must be identical in order for a consensus
                                to be called at that position. [0.7]
          --Ncutoff NCUTOFF     Maximum fraction of Ns allowed in a consensus [1.0]
          --readlength READ_LENGTH
                                Length of the input read that is being used. [80]
          --read_type READ_TYPE
                                Type of read. Options: dual_map: both reads map
                                properly. Doesn't consider read pairs where only one
                                read maps. mono_map: considers any read pair where one
                                read maps. [mono_map]
          --isize ISIZE         maximum distance between read pairs
          --read_out ROUT       How often you want to be told what the program is
                                doing. [1000000]

    DuplexMaker.py
        DCS Filter
        Version 2.0
        By Brendan Kohrn and Scott Kennedy(1)
        (1) Department of Pathology, University of Washington School of Medicine, Seattle, WA 98195 
        Based on work by Scott Kennedy, Mike Schmitt and Joe Hiatt
        October 23, 2013

        Written for Python 2.7.3
        Required modules: Pysam, Samtools, BioPython

        Inputs:
            A position-sorted paired-end BAM file containing SSCSs
            
        Outputs: 
            1: A paired-end BAM file containing DCSs
            2: A single-end BAM file containing unpaired DCSs
            3: A pair of fastq files containing DCSs for use in realligning.
            
            Note: Quality scores and cigar strings in these files are meaningless. 

        This program goes through the input file by position, making DCSs as it goes and writing them to file.  At the end of the run, any unpaired DCSs are written to a file ending in _UP.bam.  

        usage: DuplexMaker2.2.py [-h] [--infile INFILE] [--outfile OUTFILE]
                                 [--Ncutoff NCUTOFF] [--readlength READ_LENGTH]

        arguments:
          -h, --help            show this help message and exit
          --infile INFILE       input BAM file
          --outfile OUTFILE     output BAM file
          --Ncutoff NCUTOFF     Maximum percentage of Ns allowed in a consensus [1]
          --readlength READ_LENGTH
                                Length of the input read that is being used.  [80]


    count-muts.py
        Count Muts
        
        by Mike Schmitt
        Version 1.4
        October 28, 2013
        Modified by Brendan Kohrn to allow n-length indels

            
        This script pulls out the mutation frequencies from a pileup file given as stdin.

        Sites with < 20x depth, and sites with clonal mutations (defined as >30% reads mutated), are excluded from analysis by default.

        Usage:

        cat seq.pileup | python count-muts.py

    count-muys-unique.py
        Count Muts Unique
        
        by Mike Schmitt
        Version 1.3
        October 28, 2013
        Modified from count-muts.py v1.1
        Edited by Brendan Kohrn to fix a problem with 0 depth where no 0 depth should be and to allow n-length indels
            
        This script pulls out the mutation frequencies from a pileup file given as stdin.

        Sites with < 20x depth, and sites with clonal mutations (defined as >30% reads mutated), are excluded from analysis by default.

        This version counts each mutation exactly once (i.e. clonal expansions are counted as a single mutation)

        Usage:

        cat seq.pileup | python count-muts-unique.py

    mut-position.py
        Mut Position

        By Mike Schmitt
        Version 1.2
        October 28, 2013
        Modified by Brendan Kohrn to allow n-length indels

        This script gives position-specific mutation frequencies from a tagcounts file given as stdin.

        The output is tab-delimited and specifies:
        chromosome number, template base, nucleotide position, depth, mutations to T, C, G, A, insertions, deletions

        Usage:

        cat seq.pileup | python mut-position.py