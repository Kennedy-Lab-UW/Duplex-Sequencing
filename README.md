Duplex Sequencing
=================
##README

Duplex Sequencing software package  
Version 2.0  
February 20, 2014  
Programs by Scott Kennedy(1), Brendan Kohrn, and Mike Schmitt(1)  
Several steps are based on prior work by Joe Hiatt  
(1) Department of Pathology, University of Washington School of Medicine, Seattle, WA 98195

1. Glossery  
    -Single Stranded Consensus Sequence (SSCS):  
        A construct created by comparing multiple reads and deciding ambiguities by simple majority.  SSCSs are created by ConsensusMaker.py.  Quality scores attached to SSCSs are meaningless. although the cigar strings still have meaning. 
    -Duplex Consensus Sequence (DCS):  
        A construct created by comparing two SSCSs.  Quality scores and cigar strings attached to DCS sequences are meaningless, though cigar strings regain meaning after reallignment.  
        
    -Duplex tag:  
        A random sequence of nucleotides that enables the identification of sequences resulting from the same starting molecule.  
        
    -Family:  
        A group of reads that shares the same tag sequence. 

    -Read:  
        A DNA sequence which has not been compressed by *ConsensusMaker.py*.  A raw read has not yet been modified by *tag_to_header.py*, while an SMI read has.  

     

2. Summary of process  
    These programs are meant to be run in order and result in the transformation of two input FASTQ files from an Illumina sequencing run into a paired-end BAM file containing the final DCS reads.  This workflow will also generate a file containing a list of every tag that is present and how many times it occured, as well as file containing SSCSs that didn't have a mate and were unable to make a DCS (extraConsensus.bam ).  

3. Dependencies  
    The following programs and packages must be installed on your computer.  

    BWA (written with V 0.6.2)  
    Samtools (written with V 0.1.17)  
    Python (written with V 2.7.3)  
    Pysam (written with V 0.7.5)  
    BioPython (written with V 1.62)  

4. Inputs:  
	read-1-raw-data.fq  
	read-2-raw-data.fq  
 
5. Usage  
    Create a folder with both your fastq files in it.

    *PE_BASH_MAKER.py* is a script that outputs a bash script that will execute, in order, all the steps in the data processing pipeline that are needed to obtain the final DCS reads.  

    Run *PE_BASH_MAKER.py*, making sure to input the correct read length (option --rlength), using the syntax shown below. Although it is recommended that all non-optional inputs be provided, the only inputs that are truely required are --ref, --r1src, and --r2src.  

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
          --min MINMEM          Minimum members reads need to form a SSCS consensus [3]  
          --max MAXMEM          Maximum members for reads allowed to form a SSCS consensus [1000]  
          --cut CUTOFF          Mimimum percent matching for base choice in SSCS consensus [0.7]  
          --Ncut NCUT           Maxumum percent N's allowed [0.3]  
          --rlength RLENGTH     Length of a single read [84]  
          --blength BLENGTH     length of the barcode sequence on a unprocessed single read. [12]  
          --slength SLENGTH     length of the spacer sequence in a unprocessed single read.  
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
                                parallel.  This step can save a significant amount of time, but requires twice as much memory.  

    Run the bash script from the command line with:  

    ```bash
    bash scriptname.sh 3>&1 1>&2 2>&3 | tee -a log.txt   
    ```

    where scriptname.sh is the newly-created bash script (.sh file), and log.txt is the desired name of your log file.  This should run the rest of the process through to an output paired-end BAM file, copying the contents of stderr to a log file for documentation and reporting purposes.  

    Note: Do not run the bash script in a folder containing pre-existing SAM files, as it will delete them.  

    It is strongly sugested that the final sorted BAM file undergo post-processing with picard-tools-1.70/AddOrReplaceReadGroups.jar and GATK/GenomeAnalysisTK.jar, before generating statistics.  

7. Data Outputs:  
    These are only valid when using the *PE_BASH_MAKER.py* script
    \* indicates a custom string representing one of the input files.  
    
    File Description                                               | File name
    -------------------------------------------------------------- | ---------------------------------
    BAM file containing position-sorted paired-end reads:          | PE.\*.\*.bam
    BAM file containing paired-end SSCSs:                          | SSCS.\*.\*.bam
    BAM file containing unpaired SSCSs:                            | SSCS.\*.\*\_UP.bam
    BAM file containing non-mapping or otherwise bad reads:        | SSCS.\*.\*\_NM.bam
    BAM file containing good reads with less common cigar scores:  | SSCS.\*.\*\_LCC.bam
    tagcounts file:                                                | PE.\*.\*.tagcounts
    Tagstats file:                                                 | PE.\*.\*.tagstats
    Fastq files containing DCSs:                                   | DCS.\*.\*.r1.fq and PE.\*.\*.r2.fq
    BAM file containing paired-end, sorted, alligned DCSs          | DCS.\*.\*.aln.sort.bam  

8. Live Outputs  

    The file Duplex-Process-Numbers.txt describes the number of reads in each file and the live outputs from each step.    

9. Details of the individual programs in order of use (Advanced Users).  

    This information is also found at the top of each program.  

    PE_BASH_MAKER.py  
        PE Bash Maker V 1.0  
        by Brendan Kohrn  
        
        Write a bash script to run the process.    
        
        This program makes a shell script so that the user will not need to 
        enter the commands for all the programs manually.  When using it, 
        navigate to the folder with your data with all the programs in a 
        different folder.  

        usage: PE_BASH_MAKER.py [-h] [--ref REF] [--r1src R1SRC] [--r2src R2SRC]
                                [--min MINMEM] [--max MAXMEM] [--cut CUTOFF]
                                [--Ncut NCUT] [--rlength RLENGTH] [--blength BLENGTH]
                                [--slength SLENGTH] [--progInd PROGIND]
                                [--read_type READ_TYPE] [--isize ISIZE] [--absolute]
                                [--parallel]

        Arguments:
          -h, --help            show this help message and exit
          --ref REF             FASTA file containing the reference genome (genome must be properly indexed)
          --r1src R1SRC         FASTQ file containing the raw read1 data
          --r2src R2SRC         FASTQ file containing the raw read2 data
          --min MINMEM          Minimum number of reads needed to form a SSCS consensus [3]
          --max MAXMEM          Maximum number of reads that can be evaluated to make a SSCS consensus [1000]
          --cut CUTOFF          Mimimum percent matching for base choice in SSCS
                                consensus [0.7]
          --Ncut NCUT           Maxumum percent N's allowed [0.3]
          --rlength RLENGTH     Length of a single read (not overall read-pair length) after removal of barcode and spacer sequence [84]
          --blength BLENGTH     Length of the barcode sequence on a unprocessed single
                                read. [12]
          --slength SLENGTH     length of the spacer sequence in a unprocessed single
                                read.
          --progInd PROGIND     How often you want to be told what a program is doing
                                [1000000]
          --read_type READ_TYPE
                                Type of reads allowed to be considered for consensus making. d: Properly paired reads, p: Paired-end reads where
                                both reads in the pair map, but the two are not properly paired,  m: Paired-end reads where only one read in the pair maps, n: Paired-end reads where neither read maps, s: Single end mapped reads (compatibility option, not recommended for use)[dpm]
          --isize ISIZE         Optional: Maximum distance between read pairs [-1]
          --absolute            Optional: Treat the program path as an absolute path
          --parallel            Optional: Perform the alignments of both reads in
                                parallelOptional: Perform the alignments of both reads in parallel.  This is faster but requires more memory (minimum 16 GB recommended for the human genome). 


    tag_to_header.py
        Tag To Header
        Version 2.0
        By Scott Kennedy, Joe Hiatt, Brendan Kohrn and Mike Schmitt
        October 23, 2013

        Parses duplex tags from the sequence reads and appends the tag to the header (i.e. QNAME field) of each read.  Also removes the spacer region.  

        Usage: tag_to_header.py [-h] [--infile1 INFILE1] [--infile2 INFILE2]
                                [--outfile1 OUTFILE1] [--outfile2 OUTFILE2]
                                [--barcode_length BLENGTH] [--spacer_length SLENGTH]
                                [--read_out ROUT] [--adapter ADAPTERSEQ]

        Optional arguments:
          -h, --help            show this help message and exit
          --infile1 INFILE1     Name of raw FASTQ file for read 1.
          --infile2 INFILE2     Name of raw FASTQ file for read 2.
          --outfile1 OUTFILE1   Name of output file for first FASTQ read file.
          --outfile2 OUTFILE2   Name of output file for second FASTQ read file.
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
        Based on an original script by Scott Kennedy
        November 26, 2013

        Written for Python 2.7.3
        Required modules: Pysam, Samtools

        Inputs: 
            A position-sorted paired-end BAM file containing reads that have been processed by tag_to_header.py.  

        Outputs:
            1: A paired-end BAM file containing SSCSs
            2: A single-end BAM file containing unpaired SSCSs (if --read_type is 'd')
            3: A single-end BAM file containing reads with less common cigar strings
            4: A single-end BAM file containing reads not in --read_type
            5: A tagcounts file
            6: A tagstats file
            
            Note that quality scores in outputs 1, 2, and 3 are just space fillers and do not signify anything about the quality of the sequence.   

        Usage: ConsensusMaker.py [-h] [--infile INFILE] [--tagfile TAGFILE]
                                 [--outfile OUTFILE] [--rep_filt REP_FILT]
                                 [--minmem MINMEM] [--maxmem MAXMEM] [--cutoff CUTOFF]
                                 [--Ncutoff NCUTOFF] [--readlength READ_LENGTH]
                                 [--read_type READ_TYPE] [--isize ISIZE]
                                 [--read_out ROUT]

        Optional arguments:
          -h, --help            show this help message and exit
          --infile INFILE       Name of input BAM file
          --tagfile TAGFILE     Name of output tagcounts file that contains statistics on tag family size.  Used to for quality control purposes.
          --outfile OUTFILE     Name of output BAM file
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
          --readlength READ_LENGTH Length of the input read that is being used. [84]
                                
          --read_type READ_TYPE Type of reads allowed to be considered for consensus making. d: Properly paired reads, p: Paired-end reads where
                                both reads in the pair map, but the two are not properly paired,  m: Paired-end reads where only one read in the pair maps, n: Paired-end reads where neither read maps, s: Single end mapped reads (compatibility option, not recommended for use)[dpm]
          --filt                Sets which filters are used. o: Overlap filter, s: Softclipping filter, n: N filter [osn]
          --Ncutoff             with '--filt n' enabled, sets the maximum percentage of Ns allowed in a SSCS [0.3]
          --isize ISIZE         maximum distance between read pairs
          --read_out ROUT       How often you want to be told what the program is
                                doing. [1000000]

        Details of different arguments:
            --minmem and --maxmem set the range of family sizes that can be used to make a consensus sequence.  Examples use --minmem of 3 and --maxmem of 1000
                Example 1: 
                    10 reads (read length = 80) have the same tag sequence.  Of these 10, 9 of them have a CIGAR string of 80M, while one has a cigar string of 39M1I40M.  Only the 9 with a CIGAR string of 80M are sent on to be made into a SSCS.  
                Example 2:
                    3 reads (read length = 80) have the same tag sequence.  Of these, 2 have a CIGAR string of 80M, and one has a cigar string of 20M1D60M.  No SSCS results.
                Example 3: 
                    A family with over 1000 members exists.  A random sample of 1000 reads from that family is used to make a SSCS.

            --cutoff sets the strictness of the consensus making.    
                Example (--cutoff = 0.7):
                    Four reads (readlength = 10) are as follows:
                        Read 1: ACTGATACTT
                        Read 2: ACTGAAACCT
                        Read 3: ACTGATACCT
                        Read 4: ACTGATACTT
                    The resulting SSCS is:
                        ACTGATACNT

            --Ncutoff, with --filt n enabled, sets the maximum percentage of Ns allowed in a SSCS.  
                Example (--Ncutoff = .1, --readlength = 20):
                    Two SSCSs are generated as follows:
                        SSCS 1: ACGTGANCTAGTNCTNTACC
                        SSCS 2: GATCTAGTNCATGACCGATA
                    SSCS 2 passes the n filter (10%) with 1/20 = 5% Ns, while SSCS 1 does not with 3/20 = 15% Ns.

            --readlength sets the length of the reads.  If this value is set incorrectly, the program will often crash with an error message about sequence length not matching quality score length or will output an empty SSCS bam file.  

            --read_type sets which reads are considered for consensus making.  Options are: 
                d:  Paired-end reads where both reads in the pair map, and where the two are properly paired (read 2 maps in the opposite direction and on the opposite strand from read 1).  Flags are 99, 83, 163, and 147  .
                p: Paired-end reads where both reads in the pair map, but the two are not properly paired.  Flags are 97, 81, 161, 145, 129, 65, 177, and 113.
                m: Paired-end reads where only one read in the pair maps.  Flags are 181, 117, 137, 133, 73, 89, 69, and 153.
                n: Paired-end reads where neither read in the pair maps and single end unmapped reads.  Flags are 141, 77, and 4.  
                s: Single end mapped reads.  Flags are 0 and 16.
                Importantly, more than 1 option can be invoked simultaneously 

            --filt sets which filters are used.  Options are: 
                o: Overlap filter. Filters out any read pairs which overlap.  Only works on  reads of type d (see above).
                s: Softclipping filter.  Filters out any reads which have been soft-clipped in alignment.  This avoids later problems with hard-clipping.  
                n: N filter. Filters out consensus sequences with a higher percentage of Ns than the threshold imposed by --Ncutoff.  Without this option, --Ncutoff doesn't do anything.  

            --isize
                If not -1, sets the maximum distance between read 1 and read 2 for the two to not be considered unpaired.  Only works if --read_type is 'd'

    DuplexMaker.py
        DCS Filter
        Version 2.0
        By Brendan Kohrn and Scott Kennedy(1)
        (1) Department of Pathology, University of Washington School of Medicine, Seattle, WA 98195 
        Based on work by Scott Kennedy, Mike Schmitt
        October 23, 2013

        Written for Python 2.7.3
        Required modules: Pysam, Samtools, BioPython

        Inputs:
            A position-sorted paired-end BAM file containing SSCSs
            
        Outputs: 
            1: A paired-end BAM file containing DCSs
            2: A single-end BAM file containing unpaired DCSs
            3: A pair of fastq files containing DCSs for use in realigning.
            
            Note: Quality scores and CIGAR strings in these files are meaningless. 

        usage: DuplexMaker2.2.py [-h] [--infile INFILE] [--outfile OUTFILE]
                                 [--Ncutoff NCUTOFF] [--readlength READ_LENGTH]

        arguments:
          -h, --help            show this help message and exit
          --infile INFILE       input BAM file
          --outfile OUTFILE     output BAM file
          --Ncutoff NCUTOFF     Maximum percentage of Ns allowed in a consensus [1]
          --readlength READ_LENGTH
                                Length of the input read that is being used.  [84]


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
        Edited by Brendan Kohrn
            
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
