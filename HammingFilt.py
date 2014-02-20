#!/usr/bin/env python
'''
Hamming Distance Filter
Version 1.0
By Brendan Kohrn and Scott Kennedy(1)
(1) Department of Pathology, University of Washington School of Medicine, Seattle, WA 98195 
December 17, 2013

Eliminates derived tags by filtering on Hamming distance.  Anything with a hamming distance of 1 from any other tag is filtered out.  

usage: HammingFilt.py [-h] [--infile INFILE] [--read_out ROUT]

optional arguments:
  -h, --help       show this help message and exit
  --infile INFILE  input BAM file
  --read_out ROUT  how often should the program tell you how it is doing
'''

import sys
import pysam
import re
import distance
from argparse import ArgumentParser

def main():
    # Parameters to be input.
    parser=ArgumentParser()
    parser.add_argument("--infile", action="store", dest="infile", help="input BAM file", default='sys.stdin')
    parser.add_argument("--outfile", action="store", dest="outfile", help="output bam file for non-derived tags.", default = None)
    parser.add_argument('--read_out', action = 'store', dest = 'rOut', type = int, help = 'how often should the program tell you how it is doing', default=10000)
    o = parser.parse_args()
    
    if o.outfile == None:
        o.outfile = o.infile.replace('.bam', '.clean.bam')
        hamFile = o.infile.replace('.bam', '_DT.bam')
    else:
        hamFile = o.outfile.replace('.bam', '_DT.bam')

    # Initialization of all global variables, main input/output files, and main iterator and dictionaries.

    inBam = pysam.Samfile(o.infile, "rb") # Open the input BAM file
    outBam = pysam.Samfile(o.outfile, 'wb', template = inBam)
    hammingBam = pysam.Samfile(hamFile, 'wb', template = inBam)
    
    readNum = 0
    dT = 0
    hamm = 0

    fileDone=False # Initialize end of file bool
    finished=False
    readOne=True

    bamEntry = inBam.fetch( until_eof = True ) # Initialize the iterator
    firstRead = bamEntry.next() # Get the first read
    readDict = {} # Initialize the read dictionary
    firstTag=firstRead.qname.split(":")[0]

    # Start going through the input BAM file, one position at a time.

    for line in bamEntry:
        # Reinitialize first line
        readNum += 1
        if readOne==True:
            readDict[firstTag] = firstRead
            readOne=False
        
        while line.pos == firstRead.pos and fileDone==False:
            tag = line.qname.split(":")[0] #extract the barcode
            # Add the sequence to the read dictionary
            readDict[tag] = line
            try: # Keep StopIteration error from happening
                line = bamEntry.next() # Itterate the line
                readNum += 1
            except:
                fileDone = True # Tell the program that it has reached the end of the fil
                readNum += 1

            if readNum % o.rOut == 0:
                sys.stderr.write("%s reads processed\n" % (readNum))

        # Compute Hamming Distances
        FilterDict = {}
        readDictKeys=[]
        for elmt in readDict.keys():
            FilterDict[elmt] = 0
            readDictKeys.append(elmt)
        if line.pos != -1:
            for elmt in xrange(len(readDictKeys)):
                for elmt2 in xrange(elmt + 1, len(readDictKeys)):
                    HammingDist = distance.fast_comp(readDictKeys[elmt], readDictKeys[elmt2])
                    if HammingDist == 1:
                        FilterDict[readDictKeys[elmt]] = 1
                        FilterDict[readDictKeys[elmt2]] = 1
        
        firstRead = line # Store the present line for the next group of lines
        firstTag = firstRead.qname
        readOne=True
        dictKeys = [x * FilterDict[x] for x in readDictKeys]
        for elmt in dictKeys:
            if elmt != '':
                hammingBam.write(readDict[elmt])
                dT += 1
        for elmt in readDictKeys:
            FilterDict[elmt] = 1 - FilterDict[elmt]
        dictKeys = [x * FilterDict[x] for x in readDictKeys]
        for elmt in dictKeys:
            if elmt != '':
                outBam.write(readDict[elmt])
        readDict={} # Reset the read dictionary
    sys.stderr.write('Derivitive Families Removed: %s\n\n' % dT)
    print(hamm)


if __name__ == "__main__":
    main()
