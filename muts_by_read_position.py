#!/usr/bin/env python
# Original version by Brendan Kohrn
# Edited by Mike Schmitt 2/7/14 to add clonality filter

import matplotlib
matplotlib.use('Agg')
import pylab
import numpy
import sys
import re
from argparse import ArgumentParser

# myRead defines what information about each read needs to be stored.  
class myRead:
    def __init__(self, myStart, length):
        self.counts = numpy.zeros((3, length), int)
        self.fr = 'f' if myStart in ('R', 'M', 'U') else ('r' if myStart in ('r', 'm', 'u') else 'e')
        if self.fr == 'e':
            raise ValueError('first position is not a recognized pileup charecter')
        if myStart.upper() == 'M': # First position is a mutation
            self.counts[0, 0] += 1
        elif myStart.upper() == 'U': # First position is an N
            self.counts[2, 0] += 1
        elif myStart.upper() == 'R': # First position is a reference
            pass
        self.skipMe = False
        self.pos = 0
        self.closeMe = False
    
    def addMut(self):
        self.counts[0, self.pos] += 1
    
    def addIndel(self):
        self.counts[1, self.pos] += 1
    
    def addN(self):
        self.counts[2, self.pos] += 1
    
    def advance(self):
        self.pos += 1
    
    def close(self):
        if self.fr == 'r':
            return(numpy.fliplr(self.counts))
        else:
            return(self.counts)
# myCounts keeps track of the total counts so far and all current reads, as well as methods for manipulating them.
class myCounts:
    def __init__(self, length):
        self.counts = numpy.zeros((3, length), float)
        self.reads = []
        self.myLen = length
    
    def newRead(self, myStart):
        self.reads.append(myRead(myStart, self.myLen))
    
    def closeReads(self,readToClose = 0):
        closed = 0
        for read in xrange(len(self.reads)):
            if self.reads[read-closed].closeMe == True:
                self.counts += self.reads.pop(read-closed).close()
                closed += 1
    
    def advanceReads(self):
        for read in self.reads:
            if read.skipMe == False:
                read.advance()
            else:
                read.skipMe = False
    
    def muts(self):
        return(self.counts[0, :])
    
    def indels(self):
        return(self.counts[1, :])
    def ns(self):
        return(self.counts[2, :])
    
    def totals(self):
        if self.counts[0, :].sum() != 0:
            self.counts[0, :] /= self.counts[0, :].sum()
        if self.counts[1, :].sum() != 0:
            self.counts[1, :] /= self.counts[1, :].sum()
        if self.counts[2, :].sum() != 0:
            self.counts[2, :] /= self.counts[2, :].sum()

# Prepare a line for processing
def linePrep(line, maxCln):
    linebins = line.split()
    
    # Remove lines with mutations that have clonality exceeding the cutoff
    depth = max(int(linebins[3]) - linebins[4].count('N'),1)
    clonality = (float(max(linebins[4].count('T'),linebins[4].count('C'),linebins[4].count('G'),linebins[4].count('A'),linebins[4].count('t'),linebins[4].count('c'),linebins[4].count('g'),linebins[4].count('a'))) / float(depth))
    linebin = linebins[4]
    if float(clonality) > float(maxCln):
        linebin = re.sub(r'[tTcCgGaA]','.',linebin)
    
    #Convert start and end points
    linebin = re.sub(r'\$','E',linebin)
    linebin = re.sub(r'\^.\.','R', linebin)
    linebin = re.sub(r'\^.[AGCT]','M', linebin)
    linebin = re.sub(r'\^.N','U', linebin)
    linebin = re.sub(r'\^.,', 'r', linebin)
    linebin = re.sub(r'\^.[agct]','m', linebin)
    linebin = re.sub(r'\^.n','u', linebin)
    
    #Convert insertions
    newIns = map(int, re.findall(r'\+\d+', linebin))
    for length in newIns:
        rmStr = r'\+' + str(length) + "."*length
        linebin = re.sub(rmStr, str(length), linebin)
    #Convert deletions
    newDels = map(str, re.findall(r'-\d+', linebin))
    for length in newDels:
        length = int(length[1:])
        rmStr = r'-' + str(length) + "."*length
        linebin = re.sub(rmStr, 'D', linebin)
    linebin = linebin.replace('*', 'd')
    return(linebin)
   
def main():
    #Read in command-line arguments
    parser = ArgumentParser()
    parser.add_argument('-i',
                        '--infile', 
                        action ='store', 
                        dest = 'inFile', 
                        help = 'An input pileup file. If None, defaults to stdin. [None]', 
                        default = None
                        )
    parser.add_argument('-o', 
                        '--outfile', 
                        action = 'store', 
                        dest = 'outFile', 
                        help = 'A filename for the output file.  [None]', 
                        default = None
                        )
    parser.add_argument('-l',
                        '--rlength',
                        type=int,
                        action = 'store',
                        dest = 'rlength',
                        default = '84',
                        help = 'The length of a single read'
                        )
    parser.add_argument('-C',
                        '--max_clonality',
                        type = float,
                        action = 'store',
                        dest = 'max_clonality',
                        help = 'Maximum clonality to allow when considering a position [0.1]',
                        default = 0.1,
                        )
    o = parser.parse_args()
    
    # If an imput file is given, use it; otherwise, use stdin
    if o.inFile != None:
        f = open(o.inFile, 'r')
    else:
        f = sys.stdin

    counter = myCounts(o.rlength)
    linenum = 0
    lineskips = 0

    for line in f:
        linebin = linePrep(line, o.max_clonality)
        linenum += 1
        readNum = 0
        skips = 0
        if linenum % 10000 == 0:
            print('%s lines processed' % linenum)
        try:
            while readNum < len(linebin):
                # Check what the identity of a charecter is
                if linebin[readNum] in ('M', 'R', 'U', 'm', 'r', 'u'):
                    # Start a new read
                    counter.newRead(linebin[readNum]) 
                elif linebin[readNum] == 'E':
                    # Mark a read to be closed later
                    skips += 1
                    counter.reads[readNum - skips].closeMe = True
                elif linebin[readNum] in ('A', 'G', 'C', 'T', 'a', 'g', 'c', 't'):
                    # Count a mutation
                    counter.reads[readNum-skips].addMut()
                elif linebin[readNum] in ('1', '2', '3', '4', '5', '6', '7', '8', '9'):
                    # Count an indel
                    counter.reads[readNum-1-skips].addIndel()
                    # Scroll the read with the indel forward the length of the indel, and mark indel length charecters to be skipped
                    tst = 0
                    indelLength = ''
                    while readNum + tst < len(linebin) and linebin[readNum+tst] in ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9'):
                        indelLength += linebin[readNum + tst]
                        skips += 1
                        tst += 1
                    for x in xrange(int(indelLength)):
                        counter.reads[readNum-1-skips+tst].advance()
                elif linebin[readNum] == 'D':
                    counter.reads[readNum-skips-1].addIndel()
                    skips += 1
                elif linebin[readNum] == 'd':
                    counter.reads[readNum-skips].skipMe = True
                elif linebin[readNum] in ('N', 'n'):
                    # Count an N
                    counter.reads[readNum-skips].addN()
                elif linebin[readNum] in ('.', ','):
                    pass
                else:
                    raise ValueError("Not a valid pileup value at %s: %s" % (readNum, linebin[readNum]))
                readNum += 1
        except Exception:
            print('%s[%s]: %s' % (linenum, readNum, linebin))
            raise
        # Close any reads marked for closing
        counter.closeReads()
        # Advance all reads
        counter.advanceReads()
    
    # Generate and save the graphs.
    counter.totals()
    myX = range(1, o.rlength + 1)
    ax1 = pylab.subplot(3, 1, 1)
    pylab.plot(myX, counter.counts[0, :], 'b', linewidth = 2)
    pylab.title('Mutations by Position: %s' % o.inFile)
    pylab.ylabel("% Mutations")
    pylab.setp(ax1.get_xticklabels(), visible=False)
    ax2 = pylab.subplot(3, 1, 2)
    pylab.plot(myX, counter.counts[1, :], 'r', linewidth = 2)
    pylab.title('Indels by Position: %s' % o.inFile)
    pylab.ylabel("% Indels")
    pylab.setp(ax2.get_xticklabels(), visible=False)
    ax3 = pylab.subplot(3, 1, 3)
    pylab.plot(myX, counter.counts[2, :], 'k', linewidth = 2)
    pylab.title('Ns by Position: %s' % o.inFile)
    pylab.ylabel("% Ns")
    pylab.xlabel("Read Position")

    pylab.savefig(o.outFile)
    
    # Generate and save the data file
    outWrite = counter.counts.transpose()
    outFile = open(o.outFile + '.dat', 'w')
    for n in range(o.rlength):
        outStr = ""
        for m in range(3):
            outStr += ' %s' % outWrite[n, m]
        outFile.write(outStr + '\n')
    outFile.close()
    if o.inFile != None:
        f.close()

if __name__ == "__main__":
    main()
