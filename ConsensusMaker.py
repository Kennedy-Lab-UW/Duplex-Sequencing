#!/usr/bin/env python
'''
Consensus Maker
Version 2.0
By Brendan Kohrn and Scott Kennedy(1)
(1) Department of Pathology, University of Washington School of Medicine, Seattle, WA 98195
Based on work by Scott Kennedy
January 21, 2014

Written for Python 2.7.3
Required modules: Pysam, Samtools

Inputs: 
    A position-sorted paired-end BAM file containing reads with a duplex tag in the header.  

Outputs:
    1: A paired-end BAM file containing SSCSs
    2: A single-end BAM file containing unpaired SSCSs (if --read_type is 'd')
    3: A single-end BAM file containing reads with less common cigar strings
    4: A single-end BAM file containing reads not in --read_type
    5: A tagcounts file
    
    Note that quality scores in outputs 1, 2, and 3 are just space fillers and do not signify anything about the quality of the sequence.  

The program starts at the position of the first good read, determined by the type of read specified on startup.  It then goes through the file until it finds a new position, saving all reads as it goes.  When it finds a new position, it sends the saved reads to the consensus maker, one tag at a time, untill it runs out of tags.  Consensus sequences are saved until their mates come up, at which point both are written to the output BAM file, read 1 first.  After making consensuses with the reads from the first position, it continues on through the origional file until it finds another new position, sends those reads to the consensus maker, and so on until the end of the file.  At the end of the file, any remaining reads are sent through the consensus maker, and any unpaired consensuses are written to a file ending in _UP.bam.  

In the future, this program may be able to autodetect read length.  

usage: ConsensusMaker.py [-h] [--infile INFILE] [--tagfile TAGFILE]
                         [--outfile OUTFILE] [--rep_filt REP_FILT]
                         [--minmem MINMEM] [--maxmem MAXMEM] [--cutoff CUTOFF]
                         [--Ncutoff NCUTOFF] [--readlength READ_LENGTH]
                         [--read_type READ_TYPE] [--isize ISIZE]
                         [--read_out ROUT] [--filt FILT]

optional arguments:
  -h, --help            show this help message and exit
  --infile INFILE       input BAM file
  --tagfile TAGFILE     output tagcounts file
  --outfile OUTFILE     output BAM file
  --rep_filt REP_FILT   Remove tags with homomeric runs of nucleotides of
                        length x. [9]
  --minmem MINMEM       Minimum number of reads allowed to comprise a
                        consensus. [3]
  --maxmem MAXMEM       Maximum number of reads allowed to comprise a
                        consensus. [1000]
  --cutoff CUTOFF       Percentage of nucleotides at a given position in a
                        read that must be identical in order for a consensus
                        to be called at that position. [0.7]
  --Ncutoff NCUTOFF     With --filt 'n', maximum fraction of Ns allowed in a
                        consensus [1.0]
  --readlength READ_LENGTH
                        Length of the input read that is being used. [84]
  --read_type READ_TYPE
                        A string specifying which types of read to consider.
                        Read types: n: Neither read 1 or read 2 mapped. m:
                        Either read 1 or read 2 mapped, but not both. p: Both
                        read 1 and read 2 mapped, not a propper pair. d: Both
                        read 1 and read 2 mapped, propper pair. s: Single
                        ended reads ['dpm']
  --isize ISIZE         maximum distance between read pairs
  --read_out ROUT       How often you want to be told what the program is
                        doing. [1000000]
  --filt FILT           A string indicating which filters should be
                        implemented. Filters: s: Softclipping filter. o:
                        Overlap filter. n: N filter. ['osn']


Details of different arguments:
    --minmem and --maxmem set the range of family sizes (constrained by cigar score) that can be used to make a consensus sequence.  Examples use --minmem of 3 and --maxmem of 1000
        Example 1: 
            Ten reads (readlength = 80) have a particular barcode.  Of these ten, nine of them have a cigar string of 80M, while one has a cigar string of 39M1I40M.  Only the nine with a cigar string of 80M are sent on to be made into a SSCS.  
        Example 2:
            Three reads (readlength 80) have a particular barcode.  Of these, two have a cigar string of 80M, and one has a cigar string of 20M1D60M.  No SSCS results.
        Example 3: 
            A family with over 1000 members exists.  A random sample of 1000 reads from that family is used to make a SSCS.
    --cutoff sets the strictness of the consensus making.    
        Example (--cutoff = 0.7):
            Four reads (readlength = 10) are as follows:
                Read 1: ACTGATACTT
                Read 2: ACTGAAACCT
                Read 3: ACTGATACCT
                Read 4: ACTGATACTT
            The resulting SSCS is:
                ACTGATACNT
    --Ncutoff, with --filt n enabled, sets the maximum percentage of Ns allowed in a SSCS.  
        Example (--Ncutoff = .1, --readlength = 20):
            Two SSCSs are generated as follows:
                SSCS 1: ACGTGANCTAGTNCTNTACC
                SSCS 2: GATCTAGTNCATGACCGATA
            SSCS 2 passes the n filter (10%) with 1/20 = 5% Ns, while SSCS 1 does not with 3/20 = 15% Ns.
    --readlength sets the length of the reads imputed.  If this value is set incorrectly, the program will often crash with an error message about sequence length not matching quality score length, or will output an empty SSCS bam file.  
    --read_type sets which reads are considered to have 'good' flags.  Options are: 
		d:  Paired-end reads where both reads in the pair map, and where the two are properly paired (read 2 maps in the opposite direction and on the opposite strand from read 1).  Flags are 99, 83, 163, and 147  .
		p: Paired-end reads where both reads in the pair map, but the two are not properly paired.  Flags are 97, 81, 161, 145, 129, 65, 177, and 113.
		m: Paired-end reads where only one read in the pair maps.  Flags are 181, 117, 137, 133, 73, 89, 69, and 153.
		n: Paired-end reads where neither read in the pair maps, and single end unmapped reads.  Flags are 141, 77, and 4.  
		s: Single end mapped reads.  Flags are 0 and 16.  
    --filt sets which filters are used.  Options are: 
		o: Overlap filter. Filters out any read pairs which overlap.  Only works on  reads of type d (see above).
		s: Softclipping filter.  Filters out any reads which have been soft-clipped in alignment.  This avoids later problems with hard-clipping.  
		n: N filter. Filters out consensus sequences with a higher percentage of Ns than the threshold imposed by --Ncutoff.  Without this option, --Ncutoff doesn't do anything.  
    --isize
        If not -1, sets the maximum distance between read 1 and read 2 for the two to not be considered unpaired.  Only works if --read_type is 'd'

'''

import sys
import pysam
import re
import random
from collections import defaultdict
from argparse import ArgumentParser

def printRead(readIn):
    sys.stderr.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (readIn.qname, readIn.flag, readIn.tid, readIn.pos, readIn.mapq, readIn.cigar, readIn.mrnm, readIn.mpos, readIn.isize, readIn.seq, readIn.qual, readIn.tags))


def consensusMaker (groupedReadsList,  cutoff,  readLength) :
    '''The consensus maker uses a simple "majority rules" algorithm to qmake a consensus at each base position.  If no nucleotide majority reaches above the minimum theshold (--cutoff), the position is considered undefined and an 'N' is placed at that position in the read.'''
    nucIdentityList=[0, 0, 0, 0, 0, 0] # In the order of T, C, G, A, N, Total
    nucKeyDict = {0:'T', 1:'C', 2:'G', 3:'A', 4:'N'}
    consensusRead = ''

    for i in xrange(readLength) : # Count the types of nucleotides at a position in a read. i is the nucleotide index within a read in groupedReadsList
        for j in xrange(len(groupedReadsList)): # Do this for every read that comprises a SMI group. j is the read index within groupedReadsList
            try:
                if groupedReadsList[j][i] == 'T' :
                    nucIdentityList[0] += 1
                elif groupedReadsList[j][i] == 'C':
                    nucIdentityList[1] += 1
                elif groupedReadsList[j][i] == 'G':
                    nucIdentityList[2] += 1
                elif groupedReadsList[j][i] == 'A':
                    nucIdentityList[3] += 1
                elif groupedReadsList[j][i] == 'N':
                    nucIdentityList[4] += 1
                else:
                    nucIdentityList[4] += 1
                nucIdentityList[5] += 1
            except:
                break
        try:
            for j in [0, 1, 2, 3, 4] :
                if float(nucIdentityList[j])/float(nucIdentityList[5]) > cutoff :
                    consensusRead += nucKeyDict[j]
                    break
                elif j==4:
                    consensusRead += 'N'
        except:
            consensusRead += 'N'
        nucIdentityList=[0, 0, 0, 0, 0, 0] # Reset for the next nucleotide position
    return consensusRead

def tagStats(tagCountsFile):
    familySizeCounts=defaultdict( lambda: 0 )

    fIn = open(tagCountsFile, 'r')
    fOut = open(tagCountsFile.replace('.tagcounts', '.tagstats'), 'w')
    for line in fIn:
        familySizeCounts[int(line.strip().split()[1].split(":")[0])] += 1
    fIn.close()
    
    totals = 0
    for size in familySizeCounts.keys():
        familySizeCounts[size] *= int(size)
        totals += int(familySizeCounts[size])
    
    for size in sorted(familySizeCounts.keys()):
        fOut.write("%s\t%s\n" % (size, float(familySizeCounts[size])/float(totals)))
    
    fOut.close()
    return(True)

def main():
    #Parameters to be input.
    parser=ArgumentParser()
    parser.add_argument("--infile", action="store", dest="infile", help="input BAM file", required=True)
    parser.add_argument("--tagfile",  action="store",  dest="tagfile", help="output tagcounts file",  default='sys.stdout', required=True)
    parser.add_argument("--outfile",  action="store", dest="outfile", help="output BAM file", required=True)
    parser.add_argument("--rep_filt", action="store",  type=int, dest='rep_filt', help="Remove tags with homomeric runs of nucleotides of length x. [9]", default=9 )
    parser.add_argument('--minmem', type=int, default=3, dest='minmem', help="Minimum number of reads allowed to comprise a consensus. [3] ")
    parser.add_argument('--maxmem', type=int, default=1000, dest='maxmem', help="Maximum number of reads allowed to comprise a consensus. [1000]")
    parser.add_argument('--cutoff', type=float, default=.7, dest='cutoff', help="Percentage of nucleotides at a given position in a read that must be identical in order for a consensus to be called at that position. [0.7]")
    parser.add_argument('--Ncutoff', type=float, default=1, dest='Ncutoff', help="With --filt 'n', maximum fraction of Ns allowed in a consensus [1.0]")
    parser.add_argument('--readlength', type=int, default=84, dest='read_length', help="Length of the input read that is being used. [80]")
    parser.add_argument('--read_type', type=str, action="store", dest='read_type', default="dpm", help="A string specifying which types of read to consider.  Read types: n: Neither read 1 or read 2 mapped.  m: Either read 1 or read 2 mapped, but not both.  p: Both read 1 and read 2 mapped, not a propper pair.  d: Both read 1 and read 2 mapped, propper pair.  s: Single ended reads\n\t\t['dpm']")
    parser.add_argument('--isize', type = int, default=-1, dest='isize', help="maximum distance between read pairs")
    parser.add_argument('--read_out', type = int, default = 1000000, dest = 'rOut', help = 'How often you want to be told what the program is doing. [1000000]')
    parser.add_argument('--filt', type=str, default='osn', dest='filt', help="A string indicating which filters should be implemented.  Filters: s: Softclipping filter.  o: Overlap filter.  n: N filter.  ['osn']")
    o = parser.parse_args()

    # Initialization of all global variables, main input/output files, and main iterator and dictionaries.
    goodFlag=[]
    if 'd' in o.read_type:
        goodFlag.extend((99, 83, 163, 147))
    if 'm' in o.read_type:
        goodFlag.extend((181, 117, 137, 133, 73, 89, 69, 153))
    if 'p' in o.read_type:
        goodFlag.extend((97, 81, 161, 145, 129, 65, 177, 113))
    if 'n' in o.read_type:
        goodFlag.extend((141, 77, 4))
    if 's' in o.read_type:
        goodFlag.extend((0, 16))


    inBam = pysam.Samfile( o.infile, "rb" ) # Open the input BAM file
    outBam = pysam.Samfile( o.outfile, "wb", template = inBam ) # Open the output BAM file
    outNC1 = pysam.Samfile( o.outfile.replace(".bam","_LCC.bam"),"wb", template = inBam )
    nonMap = pysam.Samfile( o.outfile.replace(".bam","_NM.bam"), "wb", template = inBam ) # File for reads with strange flags
    if o.read_type == 'd':
        extraBam = pysam.Samfile(o.outfile.replace(".bam","_UP.bam"), "wb", template = inBam)

    readNum = 0
    nM = 0
    bF = 0
    oL = 0
    sC = 0
    rT = 0
    nC = 0
    
    LCC = 0
    ConMade = 0
    if o.read_type == 'd':
        UP = 0

    fileDone=False # Initialize end of file bool
    finished=False
    readOne=False

    qualScore = 'J'*o.read_length # Set a dummy quality score

    bamEntry = inBam.fetch( until_eof = True ) # Initialize the iterator
    readWin = [bamEntry.next(), ''] # Get the first read
    winPos = 0

    readDict = {} # Initialize the read dictionary
    tagDict = defaultdict( lambda: 0 ) # Initialize the tag dictionary

    consensusDict={}


#Start going through the input BAM file, one position at a time.
    for line in bamEntry:
        winPos += 1
        readWin[winPos%2] = line
        # Reinitialize first line
        if readOne==True:
            winPos -= 1
        while (readWin[winPos%2].pos == readWin[(winPos-1)%2].pos and fileDone==False and readOne==False) or readOne == True:
            if readNum % o.rOut == 0:
                sys.stderr.write("Reads processed:" + str(readNum) + "\n")
            
            try:
                tag = readWin[winPos%2].qname.split('|')[1].split('/')[0] + (":1" if readWin[winPos%2].is_read1 == True else (":2" if readWin[winPos%2].is_read2 == True else ":se"))
                tagDict[tag] += 1
            except:
                print readNum
                raise
            
            # Overlap filter: filters out overlapping reads (with --filt o)
            overlap=False
            if 'o' in o.filt:
                if readWin[winPos%2].pos < readWin[winPos%2].mpos and readWin[winPos%2].mpos < readWin[winPos%2].pos + o.read_length and int(readWin[winPos%2].flag) in (83, 99, 147, 163):
                    overlap=True
                elif readWin[winPos%2].pos > readWin[winPos%2].mpos and readWin[winPos%2].pos < readWin[winPos%2].mpos + o.read_length and int(readWin[winPos%2].flag) in (83, 99, 147, 163):
                    overlap=True
                elif readWin[winPos%2].pos==readWin[winPos%2].mpos and int(readWin[winPos%2].flag) in (83, 99, 147, 163):
                    overlap=True
            readNum +=1

            # Softclip filter: filters out softclipped reads (with --filt s)
            softClip=False
            if 's' in o.filt:
                if readWin[winPos%2].cigar != None:
                    for tupple in readWin[winPos%2].cigar:
                        if tupple[0]==4:
                            softClip=True

            # Check if the given read is good data
            if int( readWin[winPos%2].flag ) in goodFlag and overlap==False and softClip==False: 
                if ('A'*o.rep_filt in tag) or ('C'*o.rep_filt in tag) or ('G'*o.rep_filt in tag) or ('T'*o.rep_filt in tag): 
                    # Check for bad barcodes
                    nM += 1
                    nonMap.write(readWin[winPos%2])
                    rT += 1
                else :
                    # Add the sequence to the read dictionary
                    if tag not in readDict:
                        readDict[tag] = [readWin[winPos%2].flag, readWin[winPos%2].rname, readWin[winPos%2].pos, readWin[winPos%2].mrnm, readWin[winPos%2].mpos, readWin[winPos%2].isize,{str(readWin[winPos%2].cigar):[0,readWin[winPos%2].cigar]}]

                    if str(readWin[winPos%2].cigar) not in readDict[tag][6]:
                        readDict[tag][6][str(readWin[winPos%2].cigar)]=[0,readWin[winPos%2].cigar]
                    
                    readDict[tag][6][str(readWin[winPos%2].cigar)].append(readWin[winPos%2].seq)
                    readDict[tag][6][str(readWin[winPos%2].cigar)][0]+=1
            else:
                nM += 1
                nonMap.write(readWin[winPos%2])
                if int(readWin[winPos%2].flag) not in goodFlag:
                    bF += 1
                elif overlap == True:
                    oL += 1
                elif softClip == True:
                    sC += 1
            
            winPos += 1
            if readOne == False:
                try: # Keep StopIteration error from happening at the end of a file
                    readWin[winPos%2] = bamEntry.next() # Iterate the line
                except:
                    fileDone = True # Tell the program that it has reached the end of the file
            else:
                readOne = False
        else:

            # Send reads to consensusMaker
            readOne=True
            for dictTag in readDict.keys(): # Extract sequences to send to the consensus maker
                # Cigar string filtering
                cigComp={}
                for cigStr in readDict[dictTag][6].keys(): # Determin the most common cigar string
                    cigComp[cigStr]=readDict[dictTag][6][cigStr][0]
                maxCig=max(cigComp)
                if cigComp[maxCig] >= o.minmem:
                    if cigComp[maxCig] <= o.maxmem:
                        ConMade += 1
                        consensus = consensusMaker( readDict[dictTag][6][maxCig][2:],  o.cutoff,  o.read_length )
                    else:
                        ConMade += 1
                        consensus = consensusMaker(random.sample(readDict[dictTag][6][maxCig][2:], o.maxmem), o.cutoff, o.read_length)
                    
                    for cigStr in readDict[dictTag][6].keys():
                        if cigStr != maxCig:
                            for n in xrange(2, len(readDict[dictTag][6][cigStr][2:])):
                                a = pysam.AlignedRead()
                                a.qname = dictTag
                                a.flag = readDict[dictTag][0]
                                a.seq = readDict[dictTag][6][cigStr][n]
                                a.rname = readDict[dictTag][1]
                                a.pos = readDict[dictTag][2]
                                a.mapq = 255
                                a.cigar = readDict[dictTag][6][cigStr][1]
                                a.mrnm = readDict[dictTag][3]
                                a.mpos=readDict[dictTag][4]
                                a.isize = readDict[dictTag][5]
                                a.qual = qualScore  
                                outNC1.write(a)
                                LCC += 1

                    cigComp={}

                    # Filter out consensuses with too many Ns in them
                    if (consensus.count("N" )/ len(consensus) <= o.Ncutoff and 'n' in o.filt) or ('n' not in o.filt):
                        # Write a line to the consensusDictionary
                        a = pysam.AlignedRead()
                        a.qname = dictTag
                        a.flag = readDict[dictTag][0]
                        a.seq = consensus
                        a.rname = readDict[dictTag][1]
                        a.pos = readDict[dictTag][2]
                        a.mapq = 255
                        a.cigar = readDict[dictTag][6][maxCig][1]
                        a.mrnm = readDict[dictTag][3]
                        a.mpos=readDict[dictTag][4]
                        a.isize = readDict[dictTag][5]
                        a.qual = qualScore

                        # Write SSCSs to output BAM file in read pairs.
                        altTag=dictTag.replace(("1" if "1" in dictTag else "2"),("2" if "1" in dictTag else "1"))

                        if altTag in consensusDict:
                            if a.is_read1 == True:
                                outBam.write(a)
                                outBam.write(consensusDict.pop(altTag))
                            else:
                                outBam.write(consensusDict.pop(altTag))
                                outBam.write(a)
                        else:
                            consensusDict[dictTag]=a
                    else:
                        nC += 1
        readDict={} # Reset the read dictionary
        if o.read_type == 'd':
            if o.isize != -1:
                for consTag in consensusDict.keys():
                    if consensusDict[consTag].pos + o.isize < readWin[winPos%2].pos:
                        extraBam.write(consensusDict.pop(consTag))
                        UP += 1

    # Write unpaired SSCSs
    for consTag in consensusDict.keys():
        if o.read_type == 'd':
            extraBam.write(consensusDict.pop(consTag))
            UP += 1
        else:
            outBam.write(consensusDict.pop(consTag))

    # Close BAM files
    inBam.close()
    outBam.close()
    nonMap.close()
    outNC1.close()

    if o.read_type == 'd':
        extraBam.close()

    # Write summary statistics
    sys.stderr.write("Summary Statistics: \n")
    sys.stderr.write("Reads processed:" + str(readNum) + "\n")
    sys.stderr.write("Bad reads: %s\n" % nM)
    sys.stderr.write("\tReads with Bad Flags: %s\n" % bF)
    sys.stderr.write("\tOverlapping Reads: %s\n" % oL)
    sys.stderr.write("\tSoftclipped Reads: %s\n" %sC)
    sys.stderr.write("\tRepetitive Duplex Tag: %s\n" % rT)
    sys.stderr.write("Reads with Less Common Cigar Strings: %s\n" % LCC)
    sys.stderr.write("Consensuses Made: %s\n" % ConMade)
    #sys.stderr.write("Unpaired Consensuses: %s\n" % UP)
    sys.stderr.write("Consensuses with Too Many Ns: %s\n\n" % nC)

    # Write the tag counts file.
    tagFile = open( o.tagfile, "w" )
    tagFile.write ( "\n".join( [ "%s\t%d" % ( SMI, tagDict[SMI] ) for SMI in sorted( tagDict.keys(), key=lambda x: tagDict[x], reverse=True ) ] ))
    tagFile.close()
    tagStats(o.tagfile)

if __name__ == "__main__":
    main()
