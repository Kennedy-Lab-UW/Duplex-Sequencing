Duplex Sequencing
=================
##README

Duplex Sequencing software package  
Version 2.0  
August 11, 2014  
Programs by Scott Kennedy(1), Brendan Kohrn, and Mike Schmitt(1)  
Several steps are based on prior work by Joe Hiatt  
(1) Department of Pathology, University of Washington School of Medicine, Seattle, WA 98195

1. Glossery  
    -Single Stranded Consensus Sequence (SSCS):  
        A construct created by comparing multiple reads and deciding ambiguities by simple majority.  SSCSs are created by *ConsensusMaker.py*.  Quality scores attached to SSCSs are meaningless. although the cigar strings still have meaning. 
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

   Run *PE_BASH_MAKER.py*, making sure to input the correct read length (option --rlength), using the syntax shown below. Although it is recommended that all non-optional inputs be provided, the only inputs that are truely required are --ref, --r1src, --r2src, --rlength, and --runIdentifier.  Note that read_type s will not work with the default bash template.  If you want to write your own template, consult section 9.  

    ```  
    usage: PE_BASH_MAKER.py [-h] --ref REF --r1src R1SRC --r2src R2SRC --rlength RLENGTH  
                            --runIdentifier RUNID [--min MINMEM] [--max MAXMEM]  
                            [--cut CUTOFF] [--Ncut NCUT] [--blength BLENGTH]  
                            [--slength SLENGTH] [--progInd PROGIND]  
                            [--read_type READ_TYPE] [--isize ISIZE] [--filt FILT]  
                            [--repFilt REPFILT]  
                            [--template TEMPLATE]
                            
    optional arguments:  
      -h, --help            show this help message and exit  
      --ref REF             .FASTA file containing the reference genome  
      --r1src R1SRC         .fq file containing the raw read1 data  
      --r2src R2SRC         .fq file containing the raw read2 data  
      --min MINMEM          Minimum members for SSCS consensus  
      --max MAXMEM          Maximum members for SSCS consensus  
      --cut CUTOFF          Mimimum percent matching for base choice in SSCS  
                            consensus  
      --Ncut NCUT           Maxumum percent N's allowed  
      --rlength RLENGTH     Length of a single read  
      --blength BLENGTH     Length of the barcode sequence on a unprocessed single  
                            read.  
      --slength SLENGTH     Length of the spacer sequence in a unprocessed single  
                            read.  
      --progInd PROGIND     How often you want to be told what a program is doing  
      --read_type READ_TYPE  
                            A string specifying which types of read to consider.  
                            Read types: n: Neither read 1 or read 2 mapped. m:  
                            Either read 1 or read 2 mapped, but not both. p: Both  
                            read 1 and read 2 mapped, not a propper pair. d: Both  
                            read 1 and read 2 mapped, propper pair. s: Single  
                            ended reads.   
      --isize ISIZE         Optional: Maximum distance between read pairs   
      --filt FILT           A string indicating which filters should be  
                            implemented. Filters: s: Filter out softclipped reads.  
                            o: Filter out overlapping reads. n: Filter out reads  
                            with too many Ns.   
      --runIdentifier RUNID  
                            An identifier for this particular sample and  
                            sequencing run.  
      --repFilt REPFILT     Remove tags with homomeric runs of nucleotides of  
                            length x.  
      --template TEMPLATE   Template to use with bash maker. If not specified,  
                            defaults to bash_template.sh.  
    ```  
    
    The default parameters in the provided BASH script are:
    
      --min 3 --max 1000 --cut 0.7 --blength 12 --slength 5 --progInd 1000000 --read_type dpm  --isize -1 --filt os --repFilt 9 
      
    Run the bash script from the command line with:  

    ```bash
    bash runIdentifier.script.sh 3>&1 1>&2 2>&3 | tee -a runIdentifier.se.log.txt   
    ```

    where runIdentifier is the run identifier you fed to the bash maker.  This should run the rest of the process through to an output paired-end BAM file, copying the contents of stderr to a log file for documentation and reporting purposes.  
    
    It is strongly sugested that the final sorted BAM file undergo post-processing with picard-tools-1.70/AddOrReplaceReadGroups.jar and GATK/GenomeAnalysisTK.jar, before generating statistics.  

6. Data Outputs:  
    These are only valid when using the *PE_BASH_MAKER.py* script with the default template
    \* indicates the run identifier you gave *PE_BASH_MAKER.py*.  
    
    File Description                                               | File name
    -------------------------------------------------------------- | ---------------------------------
    BAM file containing position-sorted paired-end reads:          | \*.pe.bam
    BAM file containing paired-end SSCSs:                          | \*.sscs.bam
    BAM file containing unpaired SSCSs:                            | \*.sscs\_UP.bam
    BAM file containing non-mapping or otherwise bad reads:        | \*.sscs\_NM.bam
    BAM file containing good reads with less common cigar scores:  | \*.sscs\_LCC.bam
    tagcounts file:                                                | \*.pe.tagcounts
    Tagstats file:                                                 | \*.pe.tagstats
    Fastq files containing DCSs:                                   | DCS.\*.\*.r1.fq and PE.\*.\*.r2.fq
    BAM file containing paired-end, sorted, alligned DCSs          | DCS.\*.\*.aln.sort.bam  

7. Live Outputs  

    The file Duplex-Process-Numbers.txt describes the number of reads in each file and the live outputs from each step.    

8. Program Details (Advanced Users)
   
   Details of the individual programs can be found by running that program with the -h or --help options.   

9. Creating a Custom Template (Advanced Users)
   
   In order to work with the provided bash maker, all custom templates must contain the following lines before any commands are executed.  Feel free to change the default values; the bash maker just needs to have the variable names stay the same.  :  

   ```bash
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
		readLength=100
		barcodeLength=12
		spacerLength=5
		filtersSet='os'
		readTypes='dpm'
		repFilt=9
		readOut=1000000
		
		#NONDEFAULTS
		
		#FINAL_READ_LENGTH
		readLength=$((readLength-barcodeLength-spacerLength))
```
   
   Following this, the programs should be executed in the following order:  
   
   *tag_to_header.py*, *bwa*, *samtools sort*, *ConsensusMaker.py*, *samtools sort*, *DuplexMaker.py*, *bwa*  
   
   
