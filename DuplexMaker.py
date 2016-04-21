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
    1: A paired-end BAM file containing DCSs
    2: A single-end BAM file containing unpaired DCSs
    3: A pair of fastq files containing DCSs for use in realigning.
    
    Note: Quality scores and cigar strings in these files are meaningless. 

This program goes through the input file by position, making DCSs as it goes and writing them to file.  At the end of the run, any unpaired DCSs are written to a file ending in _UP.bam.  

usage: DuplexMaker.py [-h] [--infile INFILE] [--outfile OUTFILE]
                      [--Ncutoff NCUTOFF] [--readlength READ_LENGTH]
                      [--barcode_length BLENGTH] [--read_out ROUT]
                      [--gzip-fqs]

optional arguments:
  -h, --help            show this help message and exit
  --infile INFILE       input BAM file
  --outfile OUTFILE     output BAM file
  --Ncutoff NCUTOFF     Maximum percentage of Ns allowed in a consensus [1.0]
  --readlength READ_LENGTH
                        Length of the input read that is being used. [84]
  --barcode_length BLENGTH
                        Length of the duplex tag sequence. Should match the value in tag_to_header.  [12]
  --read_out ROUT       How often you want to be told what the program is
                        doing. [1000000]
  --gzip-fqs            Output gzipped fastqs [False]
'''

import sys
import pysam
import re
import gzip
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from collections import defaultdict
from argparse import ArgumentParser

def printRead(readIn):
    sys.stderr.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (readIn.qname, readIn.flag, readIn.tid, readIn.pos, readIn.mapq, readIn.cigar, readIn.mrnm, readIn.mpos, readIn.isize, readIn.seq, readIn.qual, readIn.tags))

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
    fn = outfile.replace('.bam','')+"."+end+".fq"
    if gzip_fastq:
        fn = fn + ".gz"
        return gzip.open(fn, 'wb')
    else:
        return open(fn, 'w')

def main():
    # Parameters to be input.
    parser=ArgumentParser()
    parser.add_argument("--infile", action="store", dest="infile", help="input BAM file", required=True)
    parser.add_argument("--outfile",  action="store", dest="outfile", help="output BAM file",  required=True)
    parser.add_argument('--Ncutoff', type=float, default=1.0, dest='Ncutoff', help="Maximum percentage of Ns allowed in a consensus [1.0]")
    parser.add_argument('--readlength', type=int, default=84, dest='read_length', help="Length of the input read that is being used. [84]")
    parser.add_argument('--barcode_length', type = int, default = 12, dest = 'blength', help = 'Length of the duplex tag sequence. Should match the value in tag_to_header.  [12]')
    parser.add_argument('--read_out', type = int, default = 1000000, dest = 'rOut', help = 'How often you want to be told what the program is doing. [1000000]')
    parser.add_argument('--gzip-fqs', action="store_true", default = False, dest = 'gzip_fastqs', help = 'Output gzipped fastqs [False]')
    o = parser.parse_args()

    # Initialization of all global variables, main input/output files, and main iterator and dictionaries.
    inBam = pysam.Samfile(o.infile, "rb") # Open the input BAM file
    outBam = pysam.Samfile(o.outfile, "wb", template = inBam) # Open the output BAM file
    fastqFile1 = fastq_open(o.outfile, o.gzip_fastqs, 'r1')
    fastqFile2 = fastq_open(o.outfile, o.gzip_fastqs, 'r2')

    readNum = 0
    duplexMade = 0
    uP = 0
    nC = 0

    fileDone=False # Initialize end of file bool
    finished=False
    readOne=True

    bamEntry = inBam.fetch( until_eof = True ) # Initialize the iterator
    firstRead = bamEntry.next() # Get the first read
    readDict = {} # Initialize the read dictionary
    firstTag=firstRead.qname.split(":")[0]
    qualScore = firstRead.qual # Set a dummy quality score
    consensusDict={}
    cigDum = firstRead.cigar #set a dummy cigar score

    # Start going through the input BAM file, one position at a time.
    for line in bamEntry:
        # Reinitialize first line
        readNum += 1
        if readOne==True:
            if firstRead.is_unmapped == False:
                readDict[firstTag] = [firstRead.flag, firstRead.rname, firstRead.pos, firstRead.mrnm, firstRead.mpos, firstRead.isize, firstRead.seq]
                readOne=False
        
        while line.pos == firstRead.pos and fileDone==False:
            tag = line.qname.split(":")[0] # Extract the barcode
            # Add the sequence to the read dictionary

            if line.is_unmapped == False:
                readDict[tag] = [line.flag, line.rname, line.pos, line.mrnm, line.mpos, line.isize, line.seq]
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
            firstTag = firstRead.qname.split(":")[0]
            readOne=True
            dictKeys = readDict.keys()
            
            for dictTag in readDict.keys(): # Extract sequences to send to the DCSmaker
                switchTag = dictTag[o.blength:]+dictTag[:o.blength]
                
                try:
                    consensus = DSCMaker( [readDict[dictTag][6], readDict[switchTag][6]],  o.read_length )
                    duplexMade += 1
                    # Filter out consensuses with too many Ns in them
                    if consensus.count("N")/ len(consensus) > o.Ncutoff:
                        nC += 1
                    else:
                        # Write a line to the consensusDictionary
                        a = pysam.AlignedRead()
                        a.qname = dictTag
                        a.flag = readDict[dictTag][0]
                        
                        if a.is_reverse == True:
                            tmpSeq=Seq(consensus,IUPAC.unambiguous_dna)
                            a.seq=str(tmpSeq.reverse_complement())
                        else:
                            a.seq = consensus
                        
                        a.rname = readDict[dictTag][1]
                        a.pos = readDict[dictTag][2]
                        a.mapq = 255
                        a.cigar = cigDum
                        a.mrnm = readDict[dictTag][3]
                        a.mpos=readDict[dictTag][4]
                        a.isize = readDict[dictTag][5]
                        a.qual = qualScore
                        consTag = None 
                        if dictTag in consensusDict:
                            consTag = dictTag
                        elif switchTag in consensusDict:
                            consTag = switchTag  
                        if consTag != None:
                            if a.is_read1 == True:
                                fastqFile1.write('@:%s\n%s\n+\n%s\n' %(a.qname, a.seq, a.qual))
                                outBam.write(a)
                                fastqFile2.write('@:%s\n%s\n+\n%s\n' %(consensusDict[consTag].qname, consensusDict[consTag].seq, consensusDict[consTag].qual))
                                outBam.write(consensusDict.pop(consTag))
                            else:
                                fastqFile1.write('@:%s\n%s\n+\n%s\n' %(consensusDict[consTag].qname, consensusDict[consTag].seq, consensusDict[consTag].qual))
                                outBam.write(consensusDict.pop(consTag))
                                fastqFile2.write('@:%s\n%s\n+\n%s\n' %(a.qname, a.seq, a.qual))
                                outBam.write(a)
                        else:
                            consensusDict[dictTag]=a

                    del readDict[dictTag]
                    del readDict[switchTag]
                
                except:
                    pass

        readDict={} # Reset the read dictionary

    
    # Close BAM files
    inBam.close()

    # Write unpaired DCSs
    for consTag in consensusDict.keys():
        a = pysam.AlignedRead()
        a.qname = consTag
        a.flag = 4
        a.seq = '.' * o.read_length
        a.rname = consensusDict[consTag].rname
        a.pos = consensusDict[consTag].pos
        a.mapq = 255
        a.cigar = cigDum
        a.mrnm = consensusDict[consTag].mrnm
        a.mpos=consensusDict[consTag].pos
        a.isize = consensusDict[consTag].isize
        a.qual = qualScore
        if consensusDict[consTag].is_read1 == False:
            fastqFile1.write('@:%s\n%s\n+\n%s\n' %(a.qname, a.seq, a.qual))
            outBam.write(a)
            fastqFile2.write('@:%s\n%s\n+\n%s\n' %(consensusDict[consTag].qname, consensusDict[consTag].seq, consensusDict[consTag].qual))
            outBam.write(consensusDict.pop(consTag))
        else:
            fastqFile1.write('@:%s\n%s\n+\n%s\n' %(consensusDict[consTag].qname, consensusDict[consTag].seq, consensusDict[consTag].qual))
            outBam.write(consensusDict.pop(consTag))
            fastqFile2.write('@:%s\n%s\n+\n%s\n' %(a.qname, a.seq, a.qual))
            outBam.write(a)
        uP += 1
    fastqFile1.close()
    fastqFile2.close()
    outBam.close()

    # Write summary statistics.  Duplexes made includes unpaired duplexes    
    sys.stderr.write("Summary Statistics: \n")
    sys.stderr.write("Reads Processed: %s\n" % readNum)
    sys.stderr.write("Duplexes Made: %s\n" % duplexMade)
    sys.stderr.write("Unpaired Duplexes: %s\n" % uP)
    sys.stderr.write("N-clipped Duplexes: %s\n" % nC)

if __name__ == "__main__":
    main()
