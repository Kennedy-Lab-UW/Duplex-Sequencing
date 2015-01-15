#!/bin/python
# SRAFixer.py
# by Brendan Kohrn, January 15, 2015
# 
# This program is meant to refomat reads downladed from the SRA
# (http://www.ncbi.nlm.nih.gov/sra/), including the example data set 
# (http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?run=SRR1613972).  
# It is necessary because when the reads are downloaded from the SRA, 
# they are given a different label than they had going into it, which  
# causes an error if the reads are fed straight into tag_to_header.py.
# This program fixes that error.
#
# Example: 
# tag_to_header.py is expecting the name line of a fastq to look 
# something like 
# @HWI-7001239F_017:1:1101:1226:2127/1
# but when it comes out of SRA, it looks like 
# @SRR1613972.1.1 HWI-7001239F_017:1:1101:1226:2127 length=101
#
# This program will change the read name to be
# @SRR1613972.1.1_HWI-7001239F_017:1:1101:1226:2127_length=101/1

import sys
from argparse import ArgumentParser



parser = ArgumentParser()
parser.add_argument('--infile', dest = 'infile', 
                    help = 'Input raw fastq file.  ', 
                    required=True)
parser.add_argument('--outfile', dest = 'outfile', 
                    help = 'Output fastq file.  ', 
                    required=True)
o = parser.parse_args()

infile = open(o.infile, 'r')
outfile = open(o.outfile, 'w')
readsProcessed = 0

for line in infile:
    if '@' in line or '+' in line:
        readNum = line.split(' ')[0].split('.')[2]
        repLine = "%s/%s\n" %(line.strip().replace(' ', '_'),readNum)
        outfile.write(repLine)
        if '@' in line: readsProcessed += 1
        if readsProcessed % 1000 == 0:
            print("Reads Processed: %s" % readsProcessed)
    else:
        outfile.write(line)
    
infile.close()
outfile.close()
