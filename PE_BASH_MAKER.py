#write a bash script to run the process
#This program makes a shell script so that the user will not need to enter the commands for all the programs himself.  When using it, navigate to the folder with your data, with all the programs in a different folder, and enter a relative path to that folder from the folder with your data.  This way, the program will be able to auto-build the path to the programs.  

import sys, os
import re
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("--ref", action="store", dest="ref", help=".FASTA file containing the reference genome", default='sys.stdin')
parser.add_argument("--r1src", action="store", dest='r1src', help=".fq file containing the raw read1 data", default="seq1.fq")
parser.add_argument("--r2src", action="store", dest="r2src", help=".fq file containing the raw read2 data", default="seq2.fq")
parser.add_argument("--min", action="store", dest="minMem", help="Minimum members for SSCS consensus", default="3")
parser.add_argument("--max", action="store", dest="maxMem", help="Maximum members for SSCS consensus", default="1000")
parser.add_argument("--cut", action="store", dest="cutOff", help="Mimimum percent matching for base choice in SSCS consensus", default=".8")
parser.add_argument("--Ncut", action="store", dest="Ncut", help="Maxumum percent N's allowed", default=".1")
parser.add_argument("--rlength", type=int, action="store", dest="rlength", help="Length of a single read", default="80")
parser.add_argument("--read_type", type=str, action="store", dest="read_type", default="dual_map", help="Type of read.  Options: dual_map: both reads map propperly.  Doesn't consider read pairs where only one read maps.  mono_map: considers any read pair where one read maps.")
o = parser.parse_args()

spath="".join((repr(os.getcwd()).replace("'",'')+"/"+repr(sys.argv[0]).replace("'","")))
spath=spath.replace("PE_BASH_MAKER.py","")

#Set up basic file naming protocols
r1 = o.r1src.replace(".fq", "")
r2 = o.r2src.replace(".fq", "")

#Creat output BASH file
outBash=open("PE_DCS_CALC."+r1+"."+r2+".sh", "w")
outBash.write("#!/bin/bash \n\n")


outBash.write("#Filter for reads with a properly located duplex tag, then move the tag into the header \n\n")

out1 = "head." + r1 + ".fq"
out2 = "head." + r2 + ".fq"

readlength=0

else:
	outBash.write("python " + spath + "tag_to_header.py --infile1 " + o.r1src + " --infile2 " + o.r2src + " --outfile1 " + out1 + " --outfile2 " + out2 + "\n\n")
	readlength=str(o.rlength)

aln1 = r1 + ".aln"
aln2 = r2 + ".aln"

PE = "PE." + r1 + "." + r2 + ".sam"

outBash.write("#Create a paired end file \n\n")

outBash.write("bwa aln " + o.ref + " " + out1 + " > " + aln1 + "\n")
outBash.write("bwa aln " + o.ref + " " + out2 + " > " + aln2 + "\n")
outBash.write("bwa sampe " + o.ref + " " + aln1 + " " + aln2 + " " + out1 + " " + out2 + " > " + PE + "\n\n")

PEsort = "PE." + r1 + "." + r2 + ".sort"

outBash.write("#Sort the paired end file \n\n")

outBash.write("samtools view -Sbu " + PE + " |samtools sort - " + PEsort + "\n\n")

tagF = "PE." + r1 + "." + r2 + ".tagcounts"
SSCSout = "SSCS." + r1 + "." + r2 + ".bam"

outBash.write("#Find SSCSs\n\n")

outBash.write("python " + spath + "ConsensusMaker2.2.py --infile " + PEsort + ".bam --tagfile " + tagF + " --outfile " + SSCSout + " --minmem " + o.minMem + " --maxmem " + o.maxMem + " --cutoff " + o.cutOff + " --Ncutoff " + o.Ncut + " --readlength " + readlength + " --read_type " + o.read_type + "\n\n")


SSCSsort = SSCSout.replace(".bam",".sort")

outBash.write("#Resort the SSCSs \n\n")

outBash.write("samtools view -bu " + SSCSout + " | samtools sort - " + SSCSsort + "\n\n")


DCSout = SSCSout.replace("SSCS","DCS")
outBash.write("#Find the DCSs \n\n")

outBash.write("python " + spath + "DuplexMaker2.2.py --infile " + SSCSsort + ".bam --outfile " + DCSout + " --Ncutoff " + o.Ncut + " --readlength " + readlength + "\n\n")


DCSr1 = DCSout.replace(".bam",".r1.fq")
DCSr1aln = DCSr1.replace(".fq",".aln")
DCSr2 = DCSout.replace(".bam",".r2.fq")
DCSr2aln=DCSr2.replace(".fq",".aln")
DCSaln=DCSout.replace(".bam",".aln.sam")

outBash.write("bwa aln " + o.ref + " " + DCSr1  + " > " + DCSr1aln + "\n")
outBash.write("bwa aln " + o.ref + " " + DCSr2  + " > " + DCSr2aln + "\n")
outBash.write("bwa sampe " + o.ref + " " + DCSr1aln + " " + DCSr2aln + " " + DCSr1 + " " + DCSr2 + " > " + DCSaln + "\n\n")

DCSsort = "DCS." + r1 + "." + r2 + ".aln.sort"

outBash.write("samtools view -Sbu " + DCSaln + " | samtools sort - " + DCSsort + "\n")
outBash.write("samtools index " + DCSsort + "\n")
outBash.write("samtools view -F4 " + DCSsort +  ".bam \n\n")


outBash.write("rm *.sam \n")

outBash.close()
