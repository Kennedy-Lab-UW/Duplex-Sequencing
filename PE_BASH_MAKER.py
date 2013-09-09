#write a bash script to run the process
#bash PE_DCS_Calc.txt (ref) (r1src) (r2src) (min) (max) (cut) (Ncut) (rlength)

import sys
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
parser.add_argument("--rlength", action="store", dest="rlength", help="Length of a single read", default="80")
o = parser.parse_args()
 
#Set up basic file naming protocols
r1 = o.r1src.replace(".fq", "")
r2 = o.r2src.replace(".fq", "")

#Creat output BASH file
outBash=open("PE_DCS_CALC."+r1+"."+r2+".sh", "w")
outBash.write("#!/bin/bash \n\n")


outBash.write("#Filter for reads with a properly located duplex tag, then move the tag into the header \n\n")

out1 = r1 + ".fq.smi"
out2 = r2 + ".fq.smi"

outBash.write("python tag_to_header.py --infile1 " + o.r1src + " --infile2 " + o.r2src + " --outfile1 " + out1 + " --outfile2 " + out2 + "\n\n")

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

outBash.write("python ConsensusMaker2.py --infile " + PEsort + " --tagfile " + tagF + " --outfile " + SSCSout + " --minmem " + o.minMem + " --maxmem " + o.maxMem + " --cutoff " + o.cutOff + " --Ncutoff " + o.Ncut + " --readlength " + o.rlength + "\n\n")


SSCSsort = "SSCS." + r1 + "." + r2 + ".sort"

outBash.write("#Resort the SSCSs \n\n")

outBash.write("samtools view -bu " + SSCSout + " | samtools sort - " + SSCSsort + "\n\n")


DCSout = "DCS." + r1 + "." + r2 + ".bam"

outBash.write("#Find the DCSs \n\n")

outBash.write("python DuplexMaker.py --infile " + SSCSsort + " --outfile " + DCSout + " --Ncutoff " + o.Ncut + " --readlength " + o.rlength + "\n\n")


DCSr1 = "DCS." + r1 + "." + r2 + "r1.aln"
DCSr2 = "DCS." + r1 + "." + r2 + "r2.aln"
DCSaln="DCS." + r1 + "." + r2 + ".aln.sam"

outBash.write("bwa aln -b1 " + o.ref + " " + DCSout  + " > " + DCSr1 + "\n")
outBash.write("bwa aln -b2 " + o.ref + " " + DCSout  + " > " + DCSr2 + "\n")
outBash.write("bwa sampe" + o.ref + " " + DCSr1 + " " + DCSr2 + " " + DCSout + " " + DCSout + " > " + DCSaln + "\n\n")

DCSsort = "DCS." + r1 + "." + r2 + ".aln.sort"

outBash.write("samtools view -Sbu " + DCSaln + " | samtools sort - " + DCSsort + "\n")
outBash.write("samtools index " + DCSsort + "\n")
outBash.write("samtools view -F4 " + DCSsort + ".bam \n\n")


outBash.write("rm *.sam \n")

outBash.close()

