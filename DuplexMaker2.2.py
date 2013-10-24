'''
DCS Filter
Version 2.0
By Brendan Kohrn and Scott Kennedy(1)
(1) Department of Pathology, University of Washington School of Medicine, Seattle, WA 98195 
Based on work by Scott Kennedy, Mike Schmitt and Joe Hiatt
October 23, 2013

Written for Python 2.7.3
Required modules: Pysam, Samtools, BioPython

This program is intended to be run on a paired-end BAM file, sorted by read position, which has already been through the consensus maker.  It alligns SSCS's to their switchtag, and outputs a paired-end BAM file containing Duplex Consensus Sequences (DCS's) and a BAM file containing unpaired duplex consensus sequences.  

usage: DuplexMaker2.2.py [-h] [--infile INFILE] [--outfile OUTFILE]
                         [--Ncutoff NCUTOFF] [--readlength READ_LENGTH]

arguments:
  -h, --help            show this help message and exit
  --infile INFILE       input BAM file
  --outfile OUTFILE     output BAM file
  --Ncutoff NCUTOFF     Maximum percentage of Ns allowed in a consensus [1]
  --readlength READ_LENGTH
                        Length of the input read that is being used.  [80]

'''

import sys
import pysam
import re
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from collections import defaultdict
from argparse import ArgumentParser

def printRead(readIn):
    sys.stderr.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (readIn.qname, readIn.flag, readIn.tid, readIn.pos, readIn.mapq, readIn.cigar, readIn.mrnm, readIn.mpos, readIn.isize, readIn.seq, readIn.qual, readIn.tags))


def DSCMaker (groupedReadsList,  readLength) :
    #The consensus maker uses a simple "majority rules" algorithm to qmake a consensus at each base position.  If no nucleotide majority reaches above the minimum theshold (--cutoff), the position is considered undefined and an 'N' is placed at that position in the read.
    consensusRead = ''
    for i in xrange(readLength) :#rebuild consensus read taking into account the cutoff percentage
        if groupedReadsList[0][i]==groupedReadsList[1][i]:
            consensusRead += groupedReadsList[0][i]
        else:
            consensusRead += "N"
    
    return consensusRead


def main():
    ##########################################################################################################################
    #Parameters to be input.                                                 #
    ##########################################################################################################################
    parser=ArgumentParser()
    parser.add_argument("--infile", action="store", dest="infile", help="input BAM file", default='sys.stdin')
    parser.add_argument("--outfile",  action="store", dest="outfile", help="output BAM file",  default='sys.stdout')
    parser.add_argument('--Ncutoff', type=float, default=1.0, dest='Ncutoff', help="Maximum percentage of Ns allowed in a consensus [1.0]")
    #parser.add_argument('-p', action='store_true', dest='pipe', help="Output consensus reads to stdout"  )
    parser.add_argument('--readlength', type=int, default=80, dest='read_length', help="Length of the input read that is being used. [80]")
    parser.add_argument('--read_out', type = int, default = 1000000, dest = 'rOut')
    o = parser.parse_args()

    ##########################################################################################################################
    #DCS Maker module.                                          f            #
    ##########################################################################################################################

    ##########################################################################################################################
    #Initialization of all global variables, main input/output files, and main iterator and dictionaries.            #
    ##########################################################################################################################

    inBam = pysam.Samfile(o.infile, "rb") #open the input BAM file
    outBam = pysam.Samfile(o.outfile, "wb", template = inBam) #open the output BAM file
    fastqFile1 = open(o.outfile.replace('.bam','')+".r1.fq",'w')
    fastqFile2 = open(o.outfile.replace('.bam','')+".r2.fq",'w')
    #outStd = pysam.Samfile('-', 'wb', template = inBam ) #open the stdOut writer

    #if o.pipe==False:
    #   outStd.close()
    readNum = 0
    duplexMade = 0
    uP = 0

    fileDone=False #initialize end of file bool
    finished=False
    readOne=True

    bamEntry = inBam.fetch( until_eof = True ) #initialize the iterator
    firstRead = bamEntry.next() #get the first read
    readDict = {} #initialize the read dictionary
    firstTag=firstRead.qname.split(":")[0]
    qualScore = firstRead.qual #set a dummy quality score

    ##########################################################################################################################
    #Find the first good read to serve as a start point for analysis.                                #
    ##########################################################################################################################

    consensusDict={}

    cigDum = firstRead.cigar #set a dummy cigar score

    ##########################################################################################################################
    #Start going through the input BAM file, one position at a time.                                 #
    ##########################################################################################################################

    for line in bamEntry:
        #reinitialize first line
        readNum += 1
        if readOne==True:
            readDict[firstTag] = [firstRead.flag, firstRead.rname, firstRead.pos, firstRead.mrnm, firstRead.mpos, firstRead.isize, firstRead.seq]
            readOne=False
        
        while line.pos == firstRead.pos and fileDone==False:

            tag = line.qname.split(":")[0] #extract the barcode
            #add the sequence to the read dictionary
            if tag not in readDict:
                    readDict[tag] = [line.flag, line.rname, line.pos, line.mrnm, line.mpos, line.isize, line.seq]
            #if fileDone==False:
            try: #keep StopIteration error from happening
                line = bamEntry.next() #iterate the line
                readNum += 1
            except:
                fileDone = True #tell the program that it has reached the end of the fil
                readNum += 1
            
            if readNum % o.rOut == 0:
                sys.stderr.write("%s reads processed\n" % (readNum))
        else:

    ##########################################################################################################################
    #Send reads to DCSMaker                                                #
    ##########################################################################################################################

            firstRead = line #store the present line for the next group of lines
            firstTag = firstRead.qname
            readOne=True
            dictKeys = readDict.keys()
            
            for dictTag in readDict.keys(): #extract sequences to send to the DCSmaker
                switchtag = dictTag[12:]+dictTag[:12]
                
                try:
                    consensus = DSCMaker( [readDict[dictTag][6], readDict[switchtag][6]],  o.read_length )
                    duplexMade += 1
                    #Filter out consensuses with too many Ns in them
                    if consensus.count("N" )/ len(consensus) < o.Ncutoff:
                        #write a line to the consensusDictionary
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

        ##########################################################################################################################
        #Write SSCSs to output BAM file in read pairs.                                           #
        ##########################################################################################################################
                        if dictTag in consensusDict:
                            
                            if a.is_read1 == True:
                                #if o.pipe==True:
                                #   outStd.write(a)
                                #   outStd.write(consensusDict[switchtag])UG
                                fastqFile1.write('@:%s\n%s\n+\n%s\n' %(a.qname, a.seq, a.qual))
                                outBam.write(a)
                                fastqFile2.write('@:%s\n%s\n+\n%s\n' %(consensusDict[dictTag].qname, consensusDict[dictTag].seq, consensusDict[dictTag].qual))
                                outBam.write(consensusDict.pop(dictTag))
                            else:
                                #if o.pipe==True:
                                #        outStd.write(consensusDict[switchtag])
                                #        outStd.write(a)
                                fastqFile1.write('@:%s\n%s\n+\n%s\n' %(consensusDict[dictTag].qname, consensusDict[dictTag].seq, consensusDict[dictTag].qual))
                                outBam.write(consensusDict.pop(dictTag))
                                fastqFile2.write('@:%s\n%s\n+\n%s\n' %(a.qname, a.seq, a.qual))
                                outBam.write(a)
                        else:
                            consensusDict[dictTag]=a
                    
                    del readDict[dictTag]
                    del readDict[switchtag]
                
                except:
                    pass

        readDict={} #reset the read dictionary

    ##########################################################################################################################
    #Write unpaired SSCSs to extraConsensus.bam                                          #
    ##########################################################################################################################

    extraBam=pysam.Samfile(o.outfile.replace(".bam","_UP.bam"), "wb", template = inBam)
    #close BAM files
    inBam.close()
    outBam.close()
    fastqFile1.close()
    fastqFile2.close()

    for consTag in consensusDict.keys():
        extraBam.write(consensusDict.pop(consTag))
        uP += 1
    extraBam.close()
    #outStd.close()
    
    sys.stderr.write("Summary Statistics: \n")
    sys.stderr.write("Reads Processed: %s\n" % readNum)
    sys.stderr.write("Duplexes Made: %s\n" % duplexMade)
    sys.stderr.write("Unpaired Duplexes: %s\n\n" % uP)

if __name__ == "__main__":
    main()
