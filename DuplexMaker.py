#!/usr/bin/env python
'''
DCS Filter
Version 2.0
By Brendan Kohrn, Scott Kennedy(1), and Mike Schmitt(1)
(1) Department of Pathology, University of Washington School of Medicine, Seattle, WA 98195 
Based on work by Scott Kennedy, Mike Schmitt
December 17, 2013

Written for Python 2.7.3
Required modules: Pysam, Samtools, BioPython

Inputs:
    A position-sorted paired-end BAM file containing SSCSs
    
Outputs: 
    A pair of fastq files containing DCSs for use in realigning.
    
    Note: Quality scores and cigar strings in these files are meaningless. 

This program goes through the input file by position, making DCSs as it goes and writing them to file.  At the end of the run, any unpaired DCSs are written to a file ending in _UP.bam.  

usage: DuplexMaker.py [-h] [--infile INFILE] [--outfile OUTFILE]
                      [--Ncutoff NCUTOFF] [--readlength READ_LENGTH]
                      [--barcode_length BLENGTH] [--read_out ROUT]
                      [--gzip-fqs]

optional arguments:
  -h, --help            show this help message and exit
  --infile INFILE       input BAM file
  --outfile OUTFILE     output file name prefix (output FASTQS will be prefix.r1.fq and prefix.r2.fq)
  --Ncutoff NCUTOFF     Maximum percentage of Ns allowed in a consensus [1.0]
  --readlength READ_LENGTH
                        Length of the input read that is being used. [84]
  --barcode_length BLENGTH
                        Length of the duplex tag sequence. Should match the value in tag_to_header.  [12]
  --read_out ROUT       How often you want to be told what the program is
                        doing. [1000000]
  --gzip-fqs            Output gzipped fastqs (.gz will be added to the output FASTQs) [False]
'''

import sys
import pysam
import re
import gzip
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from collections import defaultdict
from argparse import ArgumentParser

def DSCMaker (groupedReadsList,  readLength) :
    '''The Duplex maker substitutes an N if the two input sequences are not identical at a position.  '''
    consensusRead = ''
    for i in xrange(readLength) :#rebuild consensus read taking into account the cutoff percentage
        if groupedReadsList[0][i]==groupedReadsList[1][i]:
            consensusRead += groupedReadsList[0][i]
        else:
            consensusRead += "N"
    
    return consensusRead

    
def fastq_open(outfile, gzip_fastq, end):
    fn = outfile+"."+end+".fq"
    if gzip_fastq:
        fn = fn + ".gz"
        return gzip.open(fn, 'wb')
    else:
        return open(fn, 'w')

def main():
    # Parameters to be input.
    parser=ArgumentParser()
    parser.add_argument("--infile", action="store", dest="infile", help="input BAM file", required=True)
    parser.add_argument("--outfile",  action="store", dest="outfile", help="output file name prefix ",  required=True)
    parser.add_argument('--Ncutoff', type=float, default=1.0, dest='Ncutoff', help="Maximum percentage of Ns allowed in a consensus [1.0]")
    parser.add_argument('--readlength', type=int, default=84, dest='read_length', help="Length of the input read that is being used. [84]")
    parser.add_argument('--barcode_length', type = int, default = 12, dest = 'blength', help = 'Length of the duplex tag sequence. Should match the value in tag_to_header.  [12]')
    parser.add_argument('--read_out', type = int, default = 1000000, dest = 'rOut', help = 'How often you want to be told what the program is doing. [1000000]')
    parser.add_argument('--gzip-fqs', action="store_true", default = False, dest = 'gzip_fastqs', help = 'Output gzipped fastqs [False]')
    o = parser.parse_args()

    # Initialization of all global variables, main input/output files, and main iterator and dictionaries.
    inBam = pysam.Samfile(o.infile, "rb") # Open the input BAM file
    fastqFile1 = fastq_open(o.outfile, o.gzip_fastqs, 'r1')
    fastqFile2 = fastq_open(o.outfile, o.gzip_fastqs, 'r2')

    readNum = 0
    outputReadNum = 1
    duplexMade = 0
    uP = 0
    nC = 0

    fileDone=False # Initialize end of file bool
    finished=False
    readOne=True

    bamEntry = inBam.fetch( until_eof = True ) # Initialize the iterator
    firstRead = bamEntry.next() # Get the first read
    readDict = {} # Initialize the read dictionary
    firstTag=firstRead.qname.split(":")[1]
    qualScore = firstRead.qual # Set a dummy quality score
    consensusDict={}

    # Start going through the input BAM file, one position at a time.
    for line in bamEntry:
        # Reinitialize first line
        readNum += 1
        if readOne==True:
            if firstRead.is_unmapped == False:
                readDict[firstTag] = firstRead.seq
                readOne=False
        
        while line.pos == firstRead.pos and fileDone==False:
            tag = line.qname.split(":")[1] # Extract the barcode
            # Add the sequence to the read dictionary

            if line.is_unmapped == False:
                readDict[tag] = line.seq
            try: # Keep StopIteration error from happening
                line = bamEntry.next() # Itterate the line
                readNum += 1
            except:
                fileDone = True # Tell the program that it has reached the end of the file
                readNum += 1
            
            if readNum % o.rOut == 0:
                sys.stderr.write("%s reads processed\n" % readNum)
        else:
            # Send reads to DCSMaker
            firstRead = line # Store the present line for the next group of lines
            firstTag = firstRead.qname
            firstTag=firstRead.qname.split(":")[1]
            readOne=True
            dictKeys = readDict.keys()
            
            for dictTag in readDict.keys(): # Extract sequences to send to the DCSmaker
                switchtag = dictTag[o.blength:]+dictTag[:o.blength]
                
                try:
                    consensus = DSCMaker( [readDict[dictTag], readDict[switchtag]],  o.read_length )
                    duplexMade += 1
                    # Filter out consensuses with too many Ns in them
                    if consensus.count("N")/ float(len(consensus)) > o.Ncutoff:
                        nC += 1
                    else:
                        # Write a line to the consensusDictionary
                        a = pysam.AlignedRead()
                        a.qname = dictTag
                        if a.is_reverse == True:
                            tmpSeq=Seq(consensus,IUPAC.unambiguous_dna)
                            a.seq=str(tmpSeq.reverse_complement())
                        else:
                            a.seq = consensus
                        a.qual = qualScore
                        
                        # Write DCSs to output FASTQ files
                        if dictTag in consensusDict:
                            line1 = '@%d:%s\n%s\n+\n%s\n' % (outputReadNum, a.qname, a.seq, a.qual)
                            line2 = '@%d:%s\n%s\n+\n%s\n' % (outputReadNum, consensusDict[dictTag].qname, consensusDict[dictTag].seq, consensusDict[dictTag].qual)
                            outputReadNum += 1
                            if a.is_read1 == True:
                                fastqFile1.write(line1)
                                fastqFile2.write(line2)
                            else:
                                fastqFile1.write(line2)
                                fastqFile2.write(line1)
                        else:
                            consensusDict[dictTag]=a

                    del readDict[dictTag]
                    del readDict[switchtag]
                
                except:
                    pass

        readDict={} # Reset the read dictionary

    
    # Close BAM files
    inBam.close()

    # Write DCSs where only one end of a pair had a consensus.
    for consTag in consensusDict.keys():
        a = pysam.AlignedRead()
        a.qname = consTag
        a.seq = '.' * o.read_length
        a.qual = qualScore
        line1 = '@%d:%s\n%s\n+\n%s\n' % (outputReadNum, a.qname, a.seq, a.qual)
        line2 = '@%d:%s\n%s\n+\n%s\n' % (outputReadNum, consensusDict[consTag].qname, consensusDict[consTag].seq, consensusDict[consTag].qual)
        outputReadNum += 1
        if consensusDict[consTag].is_read1 == False:
            fastqFile1.write(line1)
            fastqFile2.write(line2)
        else:
            fastqFile1.write(line2)
            fastqFile2.write(line1)
        uP += 1
    fastqFile1.close()
    fastqFile2.close()

    # Write summary statistics.  Duplexes made includes unpaired duplexes    
    sys.stderr.write("Summary Statistics: \n")
    sys.stderr.write("Reads Processed: %s\n" % readNum)
    sys.stderr.write("Duplexes Made: %s\n" % duplexMade)
    sys.stderr.write("Unpaired Duplexes: %s\n" % uP)
    sys.stderr.write("N-clipped Duplexes: %s\n" % nC)

if __name__ == "__main__":
    main()
