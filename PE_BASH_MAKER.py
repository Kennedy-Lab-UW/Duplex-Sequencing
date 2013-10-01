"""write a bash script to run the process
This program makes a shell script so that the user will not need to 
enter the commands for all the programs himself.  When using it, 
navigate to the folder with your data, with all the programs in a 
different folder, and enter a relative path to that folder from the 
folder with your data.  This way, the program will be able to auto-build
the path to the programs.  
"""

import sys
import os
import re
import argparse
from argparse import ArgumentParser
import time


parser = ArgumentParser()
parser.add_argument(
        "--ref", 
        action = "store", 
        dest = "ref", 
        help = ".FASTA file containing the reference genome", 
        default = 'sys.stdin'
        )
parser.add_argument(
        "--r1src", 
        action = "store", 
        dest = 'r1src', 
        help = ".fq file containing the raw read1 data", 
        default = "seq1.fq"
        )
parser.add_argument(
        "--r2src", 
        action = "store", 
        dest = "r2src", 
        help = ".fq file containing the raw read2 data", 
        default = "seq2.fq"
        )
parser.add_argument(
        "--min", 
        action = "store", 
        dest = "minMem", 
        help = "Minimum members for SSCS consensus", 
        default = "3"
        )
parser.add_argument(
        "--max", 
        action = "store", 
        dest = "maxMem", 
        help = "Maximum members for SSCS consensus", 
        default = "1000"
        )
parser.add_argument(
        "--cut", 
        action = "store", 
        dest = "cutOff", 
        help = "Mimimum percent matching for base choice in SSCS consensus", 
        default = ".8"
        )
parser.add_argument(
        "--Ncut", 
        action = "store", 
        dest = "Ncut", 
        help = "Maxumum percent N's allowed", 
        default = ".1"
        )
parser.add_argument(
        "--rlength", 
        type = int, 
        action = "store", 
        dest = "rlength", 
        help = "Length of a single read", 
        default = "80"
        )
parser.add_argument(
        "--blength", 
        type = int, 
        action = "store", 
        dest = "blength"
        help = "length of the barcode sequence on a unprocessed single read.", 
        default = "12"
        )
parser.add_argument(
        "--slength", 
        type = int, 
        action = "store", 
        dest = "slength", 
        help = "length of the spacer sequence in a unprocessed single read.",
        default = "5"
        )
parser.add_argument(
        "--progInd", 
        type = int, 
        action = "store", 
        dest = "progInd", 
        help = "how often you want to be told what a program is doing", 
        default = "1000000"
        )
parser.add_argument(
        "--read_type", 
        type = str, 
        action = "store", 
        dest = "read_type", 
        default = "dual_map", 
        help = "Type of read.  Options: \
            dual_map: both reads map properly.  Doesn't consider read pairs \
            where only one read maps.  \
            mono_map: considers any read pair where one read maps."
            )
#BUGGY...Does not work as intended, instead runs sequentially.
parser.add_argument(
        "--parallel", 
        action = "store_true", 
        dest = "parallel", 
        help = "Perform the alignments of both reads in parallel" 
        )
o = parser.parse_args()

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

arguments=sys.argv
outBash.write("#print first few lines of the log file\n")
outBash.writeF("echo bash script made: \t%s>&2\n" % (time.ctime(time.time())))
outBash.write("echo 'python ")
for arg in arguments:
    outBash.write("%s " % arg)
outBash.write("' >&2\n")
outBash.write("echo 'chmod +x PE_DCS_CALC.*.*.sh' >&2\n")
outBash.write("echo 'bash PE_DCS_CALC.%s.%s.sh 3>&1 1>&2 2>&3 | tee -a log.txt' >&2\n\n" % (r1, r2))

outBash.write("#Filter for reads with a properly located duplex tag, ")
outBash.write("then move the tag into the header \n")

out1 = r1 + ".fq.smi"
out2 = r2 + ".fq.smi"

readlength = 0

outBash.write("echo 'tag_to_header start:' >&2\n")
outBash.write("date >&2\n")
outBash.write(
        "python %stag_to_header.py --infile1 %s --infile2 %s --outfile1 %s --outfile2 %s --barcode_length %s --spacer_length %s --read_out %s\n\n" % 
        (spath o.r1src o.r2src out1 out2 o.blength o.slength o.progInd)
        )
readlength = str(o.rlength)

aln1 = r1 + ".aln"
aln2 = r2 + ".aln"

PE = "PE." + r1 + "." + r2 + ".sam"

outBash.write("#Create a paired end file \n\n")
outBash.write("echo 'BWA start:' >&2\n")
outBash.write("date >&2\n")
if o.parallel:
        outBash.write("for read_file in " + r1 + " " + r2 + "; do bwa aln " + o.ref + " ${read_file}.fq.smi > ${read_file}.aln;done\n")

else:
        outBash.write("bwa aln " + o.ref + " " + out1 + " > " + aln1 + "\n")
        outBash.write("bwa aln " + o.ref + " " + out2 + " > " + aln2 + "\n")

outBash.write(
        "bwa sampe " + o.ref + " " + aln1 + " " + aln2 + " " + out1 + 
        " " + out2 + " > " + PE + "\n\n"
        )

PEsort = "PE." + r1 + "." + r2 + ".sort"

outBash.write("#Sort the paired end file \n\n")
outBash.write("echo 'Sort 1 start:' >&2\n")
outBash.write("date >&2\n")
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
        "python " + spath + "ConsensusMaker2.2.py --infile " + 
        PEsort + ".bam --tagfile " + tagF + " --outfile " + SSCSout + 
        " --minmem " + o.minMem + " --maxmem " + o.maxMem + " --cutoff " + 
        o.cutOff + " --Ncutoff " + o.Ncut + " --readlength " + readlength + 
        " --read_type " + o.read_type + "\n\n"
        )


SSCSsort = SSCSout.replace(".bam",".sort")

outBash.write("#Resort the SSCSs \n\n")
outBash.write("echo 'Sort 2 start:' >&2\n")
outBash.write("date >&2\n")
outBash.write(
        "samtools view -bu " + SSCSout + " | samtools sort - " + 
        SSCSsort + "\n\n"
        )


DCSout = SSCSout.replace("SSCS","DCS")
outBash.write("#Find the DCSs \n\n")
outBash.write("echo 'DuplexMaker start:' >&2\n")
outBash.write("date >&2\n")
outBash.write(
        "python " + spath + "DuplexMaker2.2.py --infile " + SSCSsort + 
        ".bam --outfile " + DCSout + " --Ncutoff " + o.Ncut + 
        " --readlength " + readlength + "\n\n"
        )


DCSr1 = DCSout.replace(".bam",".r1.fq")
DCSr1aln = DCSr1.replace(".fq",".aln")
DCSr2 = DCSout.replace(".bam",".r2.fq")
DCSr2aln = DCSr2.replace(".fq",".aln")
DCSaln = DCSout.replace(".bam",".aln.sam")
outBash.write("echo 'BWA 2 start:' >&2\n")
outBash.write("date >&2\n")
outBash.write("bwa aln " + o.ref + " " + DCSr1  + " > " + DCSr1aln + "\n")
outBash.write("bwa aln " + o.ref + " " + DCSr2  + " > " + DCSr2aln + "\n")
outBash.write(
        "bwa sampe " + o.ref + " " + DCSr1aln + " " + DCSr2aln + " " + 
        DCSr1 + " " + DCSr2 + " > " + DCSaln + "\n\n"
        )

DCSsort = "DCS." + r1 + "." + r2 + ".aln.sort"
outBash.write("echo 'Samtools sort and index start:' >&2\n")
outBash.write("date >&2\n")
outBash.write(
        "samtools view -Sbu " + DCSaln + 
        " | samtools sort - " + DCSsort + "\n"
        )
outBash.write("samtools index " + DCSsort + ".bam\n")

outBash.write("rm *.sam \n")

outBash.close()
