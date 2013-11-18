"""
PE Bash Maker V 1.0
by Brendan Kohrn
October 28, 2013

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
  --slength SLENGTH     Length of the spacer sequence in a unprocessed single
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

"""

import sys
import os
import re
import argparse
from argparse import ArgumentParser
import time

def main():
    parser = ArgumentParser()
    parser.add_argument("--ref", 
            action = "store", 
            dest = "ref", 
            help = ".FASTA file containing the reference genome"
            )
    parser.add_argument("--r1src", 
            action = "store", 
            dest = 'r1src', 
            help = ".fq file containing the raw read1 data", 
            default = "seq1.fq"
            )
    parser.add_argument("--r2src", 
            action = "store", 
            dest = "r2src", 
            help = ".fq file containing the raw read2 data", 
            default = "seq2.fq"
            )
    parser.add_argument("--min", 
            action = "store", 
            dest = "minMem", 
            help = "Minimum members for SSCS consensus [3]", 
            default = "3"
            )
    parser.add_argument("--max", 
            action = "store", 
            dest = "maxMem", 
            help = "Maximum members for SSCS consensus [1000]", 
            default = "1000"
            )
    parser.add_argument("--cut", 
            action = "store", 
            dest = "cutOff", 
            help = "Mimimum percent matching for base choice in SSCS consensus [0.8]", 
            default = ".8"
            )
    parser.add_argument("--Ncut", 
            action = "store", 
            dest = "Ncut", 
            help = "Maxumum percent N's allowed [0.1]", 
            default = ".1"
            )
    parser.add_argument("--rlength", 
            type = int, 
            action = "store", 
            dest = "rlength", 
            help = "Length of a single read [101]", 
            default = "101"
            )
    parser.add_argument("--blength", 
            type = int, 
            action = "store", 
            dest = "blength", 
            help = "Length of the barcode sequence on a unprocessed single read. [12]", 
            default = "12"
            )
    parser.add_argument("--slength", 
            type = int, 
            action = "store", 
            dest = "slength", 
            help = "Length of the spacer sequence in a unprocessed single read.",
            default = "5"
            )
    parser.add_argument("--progInd", 
            type = int, 
            action = "store", 
            dest = "progInd", 
            help = "How often you want to be told what a program is doing [1000000]", 
            default = "1000000"
            )
    parser.add_argument("--read_type", 
            type = str, 
            action = "store", 
            dest = "read_type", 
            default = "dpm", 
            help = "A string specifying which types of read to consider.  Read types: \
                    n: Neither read 1 or read 2 mapped.  \
                    m: Either read 1 or read 2 mapped, but not both.  \
                    p: Both read 1 and read 2 mapped, not a propper pair.  \
                    d: Both read 1 and read 2 mapped, propper pair.  \
                    s: Single ended reads. \
                    ['dpm']"
            )
    parser.add_argument('--isize', 
            type = int, 
            default=-1, 
            dest='isize', 
            help="Optional: Maximum distance between read pairs [-1]"
            )
    parser.add_argument('--absolute', 
            action = 'store_true', 
            dest='absolute', 
            help="Optional: Treat the program path as an absolute path"
            )
    parser.add_argument("--parallel", 
            action = "store_true", 
            dest = "parallel", 
            help = "Optional: Perform the alignments of both reads in parallel.  This is faster but requires more memory (minimum 16 GB recommended). " 
            )
    parser.add_argument('--filt', 
            action="store", 
            type=str, 
            default='osnh', 
            dest='filt', 
            help="A string indicating which filters should be implemented.  Filters: \
                    s: Filter out softclipped reads.  \
                    o: Filter out overlapping reads.  \
                    n: Filter out reads with too many Ns.  \
                    h: Filter reads based on hamming distance for derived families.  \
                    ['osnh']"
            )
    o = parser.parse_args()
    
    hamming = False
    if 'h' in o.filt:
         hamming = True
         o.filt.replace('h', '')
    
    if o.absolute:
        spath = repr(sys.argv[0]).replace("'", "").replace("PE_BASH_MAKER.py", "")
    else:
        spath = "".join(
                (repr(os.getcwd()).replace("'", '') + "/" +
                repr(sys.argv[0]).replace("'", ""))
                )
        spath = spath.replace("PE_BASH_MAKER.py", "")

    #Set up basic file naming protocols
    r1 = o.r1src.replace(".fq", "")
    r2 = o.r2src.replace(".fq", "")

    #Creat output BASH file
    outBash = open("PE_DCS_CALC." + r1 + "." + r2 + ".sh", "w")
    outBash.write("#!/bin/bash \n\n")
    outBash.write("clear\n")
    arguments=sys.argv
    outBash.write("umask 000\n\n")
    outBash.write("#print first few lines of the log file\n")
    outBash.write("echo bash script made: \t%s >&2\n" % (time.ctime(time.time())))
    outBash.write("echo 'python ")
    for arg in arguments:
        outBash.write("%s " % arg)
    outBash.write("' >&2\n")
    outBash.write("echo 'bash PE_DCS_CALC.%s.%s.sh 3>&1 1>&2 2>&3 | tee -a log.txt' >&2\n\n" % (r1, r2))

    outBash.write("#Filter for reads with a properly located duplex tag, ")
    outBash.write("then move the tag into the header \n")

    out1 = r1 + ".fq.smi"
    out2 = r2 + ".fq.smi"

    readlength = 0

    outBash.write("echo 'tag_to_header start:' >&2\n")
    outBash.write("date >&2\n")
    outBash.write(
            "echo 'python %stag_to_header.py --infile1 %s --infile2 %s --outfile1 %s --outfile2 %s --barcode_length %s --spacer_length %s --read_out %s' >&2\n\n" % 
            (spath, o.r1src, o.r2src, out1, out2, o.blength, o.slength, o.progInd)
            )
    outBash.write(
            "python %stag_to_header.py --infile1 %s --infile2 %s --outfile1 %s --outfile2 %s --barcode_length %s --spacer_length %s --read_out %s\n\n" % 
            (spath, o.r1src, o.r2src, out1, out2, o.blength, o.slength, o.progInd)
            )
    readlength = o.rlength - o.blength - o.slength

    aln1 = r1 + ".aln"
    aln2 = r2 + ".aln"

    PE = "PE." + r1 + "." + r2 + ".sam"

    outBash.write("#Create a paired end file \n\n")
    outBash.write("echo 'BWA start:' >&2\n")
    outBash.write("date >&2\n")
    if o.parallel:
            outBash.write("echo 'for read_file in " + r1 + " " + r2 + "; do' >&2\n")
            outBash.write("echo 'bwa aln " + o.ref + " ${read_file}.fq.smi > ${read_file}.aln &' >&2\n")
            outBash.write("echo 'let count+=1' >&2\n")
            outBash.write("echo '[[ $((count%2)) -eq 0 ]] && wait' >&2\n")
            outBash.write("echo 'done' >&2\n\n")
            outBash.write("for read_file in " + r1 + " " + r2 + "; do\nbwa aln " + o.ref + " ${read_file}.fq.smi > ${read_file}.aln &\nlet count+=1\n[[ $((count%2)) -eq 0 ]] && wait\ndone\n\n")

    else:
            outBash.write("echo 'bwa aln " + o.ref + " " + out1 + " > " + aln1 + "' >&2\n")
            outBash.write("bwa aln " + o.ref + " " + out1 + " > " + aln1 + "\n\n")
            outBash.write("echo 'bwa aln " + o.ref + " " + out2 + " > " + aln2 + "' >&2\n")

    outBash.write(
            "echo 'bwa sampe " + o.ref + " " + aln1 + " " + aln2 + " " + out1 + 
            " " + out2 + " > " + PE + "' >&2\n"
            )
    outBash.write(
            "bwa sampe " + o.ref + " " + aln1 + " " + aln2 + " " + out1 + 
            " " + out2 + " > " + PE + "\n\n"
            )

    PEsort = "PE." + r1 + "." + r2 + ".sort"

    outBash.write("#Sort the paired end file \n\n")
    outBash.write("echo 'Sort 1 start:' >&2\n")
    outBash.write("date >&2\n")
    
    outBash.write(
            "echo 'samtools view -Sbu " + PE + " |samtools sort - " + 
            PEsort + "' >&2\n"
            )
    outBash.write(
            "samtools view -Sbu " + PE + " |samtools sort - " + 
            PEsort + "\n\n"
            )

    tagF = "PE." + r1 + "." + r2 + ".tagcounts"
    SSCSout = "SSCS." + r1 + "." + r2 + ".bam"

    outBash.write("#Find SSCSs\n\n")
    outBash.write("echo 'ConsensusMaker start:' >&2\n")
    outBash.write("date >&2\n")
    
    outBash.write(
            "echo 'python %sConsensusMaker.py --infile %s.bam --tagfile %s --outfile %s --minmem %s --maxmem %s --cutoff %s --Ncutoff %s --readlength %s --read_type %s --filt %s --isize %s --read_out %s' >&2\n" %  
            (spath, PEsort, tagF, SSCSout, o.minMem, o.maxMem, o.cutOff, o.Ncut, readlength, o.read_type, o.filt, o.isize, o.progInd)
            )
    outBash.write(
            "python %sConsensusMaker.py --infile %s.bam --tagfile %s --outfile %s --minmem %s --maxmem %s --cutoff %s --Ncutoff %s --readlength %s --read_type %s --filt %s --isize %s --read_out %s\n\n" %  
            (spath, PEsort, tagF, SSCSout, o.minMem, o.maxMem, o.cutOff, o.Ncut, readlength, o.read_type, o.filt, o.isize, o.progInd)
            )


    SSCSsort = SSCSout.replace(".bam",".sort")

    outBash.write("#Resort the SSCSs \n\n")
    outBash.write("echo 'Sort 2 start:' >&2\n")
    outBash.write("date >&2\n")
    
    outBash.write(
            "echo 'samtools view -bu " + SSCSout + " | samtools sort - " + 
            SSCSsort + "' >&2\n"
            )
    outBash.write(
            "samtools view -bu " + SSCSout + " | samtools sort - " + 
            SSCSsort + "\n\n"
            )


    DCSout = SSCSout.replace("SSCS","DCS")
    outBash.write("#Find the DCSs \n\n")
    outBash.write("echo 'DuplexMaker start:' >&2\n")
    outBash.write("date >&2\n")
    if hamming = True:
        outBash.write(
                "echo 'python %sDuplexMaker.py --infile %s.bam --outfile %s --Ncutoff %s --readlength %s --read_out %s --hamming' >&2\n" %
                (spath, SSCSsort, DCSout, o.Ncut, readlength, o.progInd))
        outBash.write(
                "python %sDuplexMaker.py --infile %s.bam --outfile %s --Ncutoff %s --readlength %s --read_out %s --hamming >&2\n" %
                (spath, SSCSsort, DCSout, o.Ncut, readlength, o.progInd))
    else:
        outBash.write(
                "echo 'python %sDuplexMaker.py --infile %s.bam --outfile %s --Ncutoff %s --readlength %s --read_out %s' >&2\n" %
                (spath, SSCSsort, DCSout, o.Ncut, readlength, o.progInd))
        outBash.write(
                "python %sDuplexMaker.py --infile %s.bam --outfile %s --Ncutoff %s --readlength %s --read_out %s >&2\n" %
                (spath, SSCSsort, DCSout, o.Ncut, readlength, o.progInd))


    DCSr1 = DCSout.replace(".bam",".r1.fq")
    DCSr1aln = DCSr1.replace(".fq",".aln")
    DCSr2 = DCSout.replace(".bam",".r2.fq")
    DCSr2aln = DCSr2.replace(".fq",".aln")
    DCSaln = DCSout.replace(".bam",".aln.sam")
    outBash.write("echo 'BWA 2 start:' >&2\n")
    outBash.write("date >&2\n")
    outBash.write("echo 'bwa aln " + o.ref + " " + DCSr1  + " > " + DCSr1aln + "' >&2\n")
    outBash.write("bwa aln " + o.ref + " " + DCSr1  + " > " + DCSr1aln + "\n\n")
    outBash.write("echo 'bwa aln " + o.ref + " " + DCSr2  + " > " + DCSr2aln + "' >&2\n")
    outBash.write("bwa aln " + o.ref + " " + DCSr2  + " > " + DCSr2aln + "\n\n")
    outBash.write(
            "echo 'bwa sampe " + o.ref + " " + DCSr1aln + " " + DCSr2aln + " " + 
            DCSr1 + " " + DCSr2 + " > " + DCSaln + "' >&2\n"
            )
    outBash.write(
            "bwa sampe " + o.ref + " " + DCSr1aln + " " + DCSr2aln + " " + 
            DCSr1 + " " + DCSr2 + " > " + DCSaln + "\n\n"
            )

    DCSsort = "DCS." + r1 + "." + r2 + ".aln.sort"
    outBash.write("echo 'Samtools sort and index start:' >&2\n")
    outBash.write("date >&2\n")
    outBash.write(
            "echo 'samtools view -Sbu " + DCSaln + 
            " | samtools sort - " + DCSsort + "' >&2\n"
            )
    outBash.write(
            "samtools view -Sbu " + DCSaln + 
            " | samtools sort - " + DCSsort + "\n\n"
            )
    outBash.write("echo 'samtools index " + DCSsort + ".bam' >&2\n")
    outBash.write("samtools index " + DCSsort + ".bam\n\n")

    outBash.write("echo 'rm *.sam' >&2\n")
    outBash.write("rm *.sam \n")

    outBash.close()

if __name__ == "__main__":
    main()
