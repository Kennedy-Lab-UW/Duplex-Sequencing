"""
PE Bash Maker V 1.11
by Brendan Kohrn
August 11, 2014

Write a bash script to run the process.  

This program makes a shell script so that the user will not need to 
enter the commands for all the programs himself.  When using it, 
navigate to the folder with your data, with all the programs in a 
different folder.  Options --ref, --r1src, --r2src, and --runIdentifier are required.  
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
            help = ".FASTA file containing the reference genome", 
            required=True
            )
    parser.add_argument("--r1src", 
            action = "store", 
            dest = "r1src", 
            help = ".fq file containing the raw read1 data", 
            required=True
            )
    parser.add_argument("--r2src", 
            action = "store", 
            dest = "r2src", 
            help = ".fq file containing the raw read2 data", 
            required=True
            )
    parser.add_argument("--min", 
            action = "store", 
            dest = "minMem", 
            help = "Minimum members for SSCS consensus"
            )
    parser.add_argument("--max", 
            action = "store", 
            dest = "maxMem", 
            help = "Maximum members for SSCS consensus"
            )
    parser.add_argument("--cut", 
            action = "store", 
            dest = "cutOff", 
            help = "Mimimum percent matching for base choice in SSCS consensus "
            )
    parser.add_argument("--Ncut", 
            action = "store", 
            dest = "Ncut", 
            help = "Maxumum percent N's allowed"
            )
    parser.add_argument("--rlength", 
            action = "store", 
            dest = "rlength", 
            help = "Length of a single read",
            required = True
            )
    parser.add_argument("--blength", 
            action = "store", 
            dest = "blength", 
            help = "Length of the barcode sequence on a unprocessed single read. "
            )
    parser.add_argument("--slength",  
            action = "store", 
            dest = "slength", 
            help = "Length of the spacer sequence in a unprocessed single read."
            )
    parser.add_argument("--progInd", 
            action = "store", 
            dest = "progInd", 
            help = "How often you want to be told what a program is doing"
            )
    parser.add_argument("--read_type", 
            type = str, 
            action = "store", 
            dest = "read_type", 
            help = "A string specifying which types of read to consider.  Read types: \
                    n: Neither read 1 or read 2 mapped.  \
                    m: Either read 1 or read 2 mapped, but not both.  \
                    p: Both read 1 and read 2 mapped, not a propper pair.  \
                    d: Both read 1 and read 2 mapped, propper pair.  \
                    s: Single ended reads. \
                    ['dpm']"
            )
    parser.add_argument("--isize",
            action = "store", 
            dest="isize", 
            help="Optional: Maximum distance between read pairs"
            )
    parser.add_argument("--filt",
            action="store", 
            dest="filt", 
            help="A string indicating which filters should be implemented.  Filters: \
                    s: Filter out softclipped reads.  \
                    o: Filter out overlapping reads.  \
                    n: Filter out reads with too many Ns.  \
                    ['osn']"
            )
    parser.add_argument("--runIdentifier",
            dest = "runID",
            help = "An identifier for this particular sample and sequencing run.", 
            required=True
            )
    parser.add_argument("--repFilt",
            action = "store",
            dest = "repFilt",
            help = "Remove tags with homomeric runs of nucleotides of length x."
            )
    parser.add_argument("--template", 
            action = "store", 
            dest="template",
            help="Template to use with bash maker.  If not specified, defaults to bash_template.sh."
            )
    o = parser.parse_args()
    outBash = open(o.runID + ".script.sh", "w")
    
    spath = repr(sys.argv[0]).replace("'", "").replace("/PE_BASH_MAKER.py", "")
    if o.template:
        inBash = open(o.template, "r")
    else:
        inBash = open(spath + "/bash_template.sh", "r")
    
    for line in inBash:
        if line.strip() != "#NONDEFAULTS":
            outBash.write(line)
        else:
            outBash.write(line + "\n")
            outBash.write("DSpath='" + spath + "'\n")
            if o.ref:
                outBash.write("alignRef='" + o.ref + "'\n")
            outBash.write("runIdentifier=" + o.runID + "\n")
            if o.r1src:
                outBash.write("read1in=" + o.r1src + "\n")
            if o.r2src:
                outBash.write("read2in=" + o.r2src + "\n")
            if o.isize:
                outBash.write("iSize=" + o.isize + "\n")
            if o.minMem:
                outBash.write("minMem=" + o.minMem + "\n")
            if o.maxMem:
                outBash.write("maxMem=" + o.maxMem + "\n")
            if o.cutOff:
                outBash.write("cutOff=" + o.cutOff + "\n")
            if o.Ncut:
                outBash.write("nCutOff=" + o.Ncut + "\n")
            if o.rlength:
                outBash.write("readLength=" + o.rlength + "\n")
            if o.blength:
                outBash.write("barcodeLength=" + o.blength + "\n")
            if o.slength:
                outBash.write("spacerLength=" + o.slength + "\n")
            if o.filt:
                outBash.write("filtersSet=" + o.filt + "\n")
            if o.read_type:
                outBash.write("readTypes=" + o.read_type + "\n")
            if o.repFilt:
                outBash.write("repFilt=" + o.repFilt + "\n")
            if o.progInd:
                outBash.write("readOut=" + o.progInd + "\n")
    
    outBash.write("\n\n#Generated on " + time.ctime(time.time()) + " with the command: \n")
    outBash.write("#python ")
    arguments=sys.argv
    for arg in arguments:
        outBash.write("%s " % arg)
    outBash.write("\n")
    outBash.close()
    inBash.close()
    
    print("Script made.  Run with:\n%s.script.sh\n" % o.runID)

if __name__ == "__main__":
    main()
