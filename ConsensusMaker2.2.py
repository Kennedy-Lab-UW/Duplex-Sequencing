'''
Consensus Maker
Version 2.0
By Brendan Kohrn and Scott Kennedy(1)
(1) Department of Pathology, University of Washington School of Medicine, Seattle, WA 98195
Based on work by Mike Schmitt and Joe Hiatt
October 23, 2013

Written for Python 2.7.3
Required modules: Pysam, Samtools

This program is intended to be run on a paired-end BAM file, sorted by read position, with duplex tags in the header and constant read length.  It will output a paired-end BAM file with single strand consensus sequences (SSCSs), and a .tagcounts file which contains the different tags (on both strands) and how many times they occur, even if they are not used in SSCS generation, in order by read.  In addition, it will output a BAM file of SSCSs which are unpaired, either because one of the pair didn't match the criteria for allignment, or because of some other reason, and a BAM file of all unconsidered sequences in the original file.  Quality scores on the output BAM files are meaningless.  The file produced by this program is meant to continue on through the duplex maker.  

The program starts at the position of the first good read, determined by the type of read specified on startup.  It then goes through the file until it finds a new position, saving all reads as it goes.  When it finds a new position, it sends the saved reads to the consensus maker, one tag at a time, untill it runs out of tags.  Consensus sequences are saved until their mates come up, at which point both are written to the output BAM file, first read first.  After emptying the reads from the first position, it continues on through the origional file until it finds another new position, sends those reads to the consensus maker, and so on until the end of the file.  At the end of the file, any remaining reads are sent through the consensus maker, and any unpaired consensuses are written to extraConsensus.bam.  

In the future, the program may be able to autodetect read length.  

usage: ConsensusMaker2.2.py [-h] [--infile INFILE] [--tagfile TAGFILE]
                            [--outfile OUTFILE] [--rep_filt REP_FILT]
                            [--minmem MINMEM] [--maxmem MAXMEM]
                            [--cutoff CUTOFF] [--Ncutoff NCUTOFF]
                            [--readlength READ_LENGTH] [--read_type READ_TYPE]

arguments:
  -h, --help            show this help message and exit
  --infile INFILE       input BAM file
  --tagfile TAGFILE     output tagcounts file
  --outfile OUTFILE     output BAM file
  --rep_filt REP_FILT   Remove tags with homomeric runs of nucleotides of
                        length x [9]
  --minmem MINMEM       Minimum number of reads allowed to comprise a
                        consensus. [0]
  --maxmem MAXMEM       Maximum number of reads allowed to comprise a
                        consensus. [100]
  --cutoff CUTOFF       Percentage of nucleotides at a given position in a
                        read that must be identical in order for a consensus
                        to be called at that position. [0]
  --Ncutoff NCUTOFF     Maximum percentage of Ns allowed in a consensus [1]
  --readlength READ_LENGTH
                        Length of the input read that is being used. [80]
  --read_type READ_TYPE
                        Type of read. 
                        Options: 
                            dual_map: both reads map propperly.  Doesn't consider read pairs where only one read maps. 
                            mono_map: considers any read pair where one read maps. 
'''

'''
ChangeLog in this version:
Inserted try in consensusMaker.  
Added cigarDictionary as readDict[6]
Added cigar string comparison before sending to consensus maker.
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
    seqDict = {}
    consensusRead = ''

    for i in xrange(readLength) : #Count the types of nucleotides at a position in a read. i is the nucleotide index within a read in groupedReadsList
        for j in xrange(len(groupedReadsList)): #Do this for every read that comprises a SMI group. j is the read index within groupedReadsList
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
                nucIdentityList[5] += 1
                seqDict[i] = nucIdentityList
            except:
                seqDict[i] = nucIdentityList
                nucIdentityList=[0, 0, 0, 0, 0, 0]
                break
        nucIdentityList=[0, 0, 0, 0, 0, 0] #reset for the next nucleotide position

    for i in xrange(readLength) :#rebuild consensus read taking into account the cutoff percentage
        try:
            for j in [0, 1, 2, 3, 4] :
                if float(seqDict[i][j])/float(seqDict[i][5]) > cutoff :
                    consensusRead += nucKeyDict[j]
                    break
                elif j==4:
                    consensusRead += 'N'
        except:
            consensusRead += 'N'
    return consensusRead

def main():
    ##########################################################################################################################
    #Parameters to be input.                                                 #
    ##########################################################################################################################
    parser=ArgumentParser()
    parser.add_argument("--infile", action="store", dest="infile", help="input BAM file", default='sys.stdin')
    parser.add_argument("--tagfile",  action="store",  dest="tagfile", help="output tagcounts file",  default='sys.stdout')
    parser.add_argument("--outfile",  action="store", dest="outfile", help="output BAM file",  default='sys.stdout')
    parser.add_argument("--rep_filt", action="store",  type=int, dest='rep_filt', help="Remove tags with homomeric runs of nucleotides of length x. [9]", default=9 )
    parser.add_argument('--minmem', type=int, default=3, dest='minmem', help="Minimum number of reads allowed to comprise a consensus. [3] ")
    parser.add_argument('--maxmem', type=int, default=1000, dest='maxmem', help="Maximum number of reads allowed to comprise a consensus. [1000]")
    parser.add_argument('--cutoff', type=float, default=.7, dest='cutoff', help="Percentage of nucleotides at a given position in a read that must be identical in order for a consensus to be called at that position. [0.7]")
    parser.add_argument('--Ncutoff', type=float, default=1, dest='Ncutoff', help="Maximum fraction of Ns allowed in a consensus [1.0]")
    parser.add_argument('--readlength', type=int, default=80, dest='read_length', help="Length of the input read that is being used. [80]")
    parser.add_argument('--read_type', action="store", dest='read_type', default="mono_map", help="Type of read.  Options: dual_map: both reads map properly.  Doesn't consider read pairs where only one read maps.  mono_map: considers any read pair where one read maps. [mono_map]")
    parser.add_argument('--isize', type = int, default=-1, dest='isize', help="maximum distance between read pairs")
    parser.add_argument('--read_out', type = int, default = 1000000, dest = 'rOut')
    o = parser.parse_args()

    ##########################################################################################################################
    #Consensus Maker module.  Takes a list of sequences, and finds a consensus by simple majority.               #
    ##########################################################################################################################


    ##########################################################################################################################
    #Initialization of all global variables, main input/output files, and main iterator and dictionaries.            #
    ##########################################################################################################################
    goodFlag=[]
    if o.read_type == "dual_map":
        goodFlag=[83, 99, 147, 163]
    elif o.read_type == "mono_map":
        goodFlag=[153, 89, 117, 181, 83, 99, 171, 163]

    inBam = pysam.Samfile( o.infile, "rb" ) #open the input BAM file
    outBam = pysam.Samfile( o.outfile, "wb", template = inBam ) #open the output BAM file
    outNC1 = pysam.Samfile( o.outfile.replace(".bam","_LCC.bam"),"wb", template = inBam )
    nonMap = pysam.Samfile( o.outfile.replace(".bam","_NM.bam"), "wb", template = inBam ) #file for reads with strange flags
    extraBam = pysam.Samfile(o.outfile.replace(".bam","_UP.bam"), "wb", template = inBam)
    #outStd = pysam.Samfile('-', 'wb', template = inBam ) #open the stdOut writer

    readNum = 0
    nM = 0
    bF = 0
    oL = 0
    sC = 0
    rT = 0
    
    LCC = 0
    ConMade = 0
    UP = 0

    fileDone=False #initialize end of file bool
    finished=False
    readOne=False

    qualScore = 'J'*o.read_length #set a dummy quality score

    bamEntry = inBam.fetch( until_eof = True ) #initialize the iterator
    readWin = [bamEntry.next(), ''] #get the first read
    winPos = 0

    readDict = {} #initialize the read dictionary
    tagDict = defaultdict( lambda: 0 ) #initialize the tag dictionary

    consensusDict={}

    ##########################################################################################################################
    #Start going through the input BAM file, one position at a time.                                 #
    ##########################################################################################################################
    for line in bamEntry:
        winPos += 1
        readWin[winPos%2] = line
        #reinitialize first line
        if readOne==True:
            winPos -= 1
        while (readWin[winPos%2].pos == readWin[(winPos-1)%2].pos and fileDone==False and readOne==False) or readOne == True:
            if readNum % o,rOut == 0:
                sys.stderr.write("Reads processed:" + str(readNum) + "\n")
            overlap=False
            if readWin[winPos%2].pos < readWin[winPos%2].mpos and readWin[winPos%2].mpos < readWin[winPos%2].pos + o.read_length and int(readWin[winPos%2].flag) in (83, 99, 147, 163):
                overlap=True
            elif readWin[winPos%2].pos > readWin[winPos%2].mpos and readWin[winPos%2].pos < readWin[winPos%2].mpos + o.read_length and int(readWin[winPos%2].flag) in (83, 99, 147, 163):
                overlap=True
            elif readWin[winPos%2].pos==readWin[winPos%2].mpos and int(readWin[winPos%2].flag) in (83, 99, 147, 163):
                overlap=True
            readNum +=1

            softClip=False
            if readWin[winPos%2].cigar != None:
                for tupple in readWin[winPos%2].cigar:
                    if tupple[0]==4:
                        softClip=True

            tag = readWin[winPos%2].qname.split('#')[1] + (":1" if readWin[winPos%2].is_read1 == True else (":2" if readWin[winPos%2].is_read2 == True else ":se"))
            tagDict[tag] += 1

            if int( readWin[winPos%2].flag ) in goodFlag and overlap==False and softClip==False: #check if the given read is good data
                if ('A'*o.rep_filt in tag) or ('C'*o.rep_filt in tag) or ('G'*o.rep_filt in tag) or ('C'*o.rep_filt in tag): 
                    #check for bad barcodes
                    nM += 1
                    nonMap.write(readWin[winPos%2])
                    rT += 1
                else :
                    #add the sequence to the read dictionary
                    if tag not in readDict:
                        readDict[tag] = [readWin[winPos%2].flag, readWin[winPos%2].rname, readWin[winPos%2].pos, readWin[winPos%2].mrnm, readWin[winPos%2].mpos, readWin[winPos%2].isize,{str(readWin[winPos%2].cigar):[0,readWin[winPos%2].cigar]}]

                    if str(readWin[winPos%2].cigar) not in readDict[tag][6]:
                        readDict[tag][6][str(readWin[winPos%2].cigar)]=[0,readWin[winPos%2].cigar]
                    
                    readDict[tag][6][str(readWin[winPos%2].cigar)].append(readWin[winPos%2].seq)
                    readDict[tag][6][str(readWin[winPos%2].cigar)][0]+=1
            else:
                nM += 1
                nonMap.write(readWin[winPos%2])
                if int(readWin[winPos%2].flag() not in goodFlag:
                    bF += 1
                elif overlap == True:
                    oL += 1
                elif softClip == True:
                    sC += 1
            
            winPos += 1
            if readOne == False:
                try: #keep StopIteration error from happening
                    readWin[winPos%2] = bamEntry.next() #iterate the line
                except:
                    fileDone = True #tell the program that it has reached the end of the file
            else:
                readOne = False
        else:

    ##########################################################################################################################
    #Send reads to consensusMaker                                                #
    ##########################################################################################################################

            readOne=True
            for dictTag in readDict.keys(): #extract sequences to send to the consensus maker
                cigComp={}
                
                for cigStr in readDict[dictTag][6].keys(): #determin the most common cigar string
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
                                a.qname = dictTag.split(':')[0]
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

                    #Filter out consensuses with too many Ns in them
                    if consensus.count("N" )/ len(consensus) < o.Ncutoff:
                        #write a line to the consensusDictionary
                        a = pysam.AlignedRead()
                        a.qname = dictTag.split(':')[0]
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

    ##########################################################################################################################
    #Write SSCSs to output BAM file in read pairs.                                           #
    ##########################################################################################################################

                        altTag=dictTag.replace(("1" if "1" in dictTag else "2"),("2" if "1" in dictTag else "1"))

                        if altTag in consensusDict:
                            if a.is_read1 == True:
                                #if o.pipe==True:
                                #   outStd.write(a)
                                #   outStd.write(consensusDict[dictTag])
                                outBam.write(a)
                                outBam.write(consensusDict.pop(altTag))
                            else:
                                #if o.pipe==True:
                                #       outStd.write(consensusDict[dictTag])
                                #       outStd.write(a)
                                outBam.write(consensusDict.pop(altTag))
                                outBam.write(a)
                        else:
                            consensusDict[dictTag]=a

        readDict={} #reset the read dictionary
        
        if o.isize != -1:
            for consTag in consensusDict.keys():
                if consensusDict[consTag].pos + o.isize < readWin[winPos%2].pos:
                    extraBam.write(consensusDict.pop(consTag))
                    UP += 1

    ##########################################################################################################################
    #Write unpaired SSCSs to extraConsensus.bam                                          #
    ##########################################################################################################################

    #close BAM files
    inBam.close()
    outBam.close()
    nonMap.close()
    outNC1.close()

    for consTag in consensusDict.keys():
        extraBam.write(consensusDict.pop(consTag))
        UP += 1

    extraBam.close()
    #outStd.close()
    ##########################################################################################################################
    #Write the tag counts file.                                                  #
    ##########################################################################################################################
    
    sys.stderr.write("Summary Statistics: \n")
    sys.stderr.write("Reads processed:" + str(readNum) + "\n")
    sys.stderr.write("Bad reads: %s\n" % nM)
    sys.stderr.write("\tNon-mapping reads: %s\n" % bF)
    sys.stderr.write("\tOverlapping Reads: %s\n" % oL)
    sys.stderr.write("\tSoftclipped Reads: %s\n" %sC)
    sys.stderr.write("\tRepetitive Duplex Tag: %s\n" % rT)
    sys.stderr.write("Reads with Less Common Cigar Strings: %s\n" % LCC)
    sys.stderr.write("Consensuses Made: %s\n" % ConMade)
    sys.stderr.write("Unpaired Consensuses: %s\n\n" % UP)

    tagFile = open( o.tagfile, "w" )
    tagFile.write ( "\n".join( [ "%s\t%d" % ( SMI, tagDict[SMI] ) for SMI in sorted( tagDict.keys(), key=lambda x: tagDict[x], reverse=True ) ] ))
    tagFile.close()

if __name__ == "__main__":
    main()
