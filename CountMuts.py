"""
Count Muts Unique

by Brendan Kohrn and Mike Schmitt
Version 1.4
October 28, 2013
Modified from count-muts.py and count-muts-unique.py
Edited by Brendan Kohrn to fix a problem with 0 depth where no 0 depth should be and to allow n-length indels, and to merge the two conditions into one program.  
This script pulls out the mutation frequencies from a pileup file given as stdin.

Sites with < 20x depth, and sites with clonal mutations (defined as >30% reads mutated), are excluded from analysis by default.

If -u is specified, this program counts each mutation exactly once (i.e. clonal expansions are counted as a single mutation)

Usage:

cat seq.pileup | python CountMuts.py
"""

from argparse import ArgumentParser
import sys
import re
from math import sqrt


def Wilson(positive,  total) :
    
        if total == 0:
            
            return 0
            
        freq = float(positive)/float(total)
        z = 1.96 #1.96 = 95%
        phat = float(positive) / total
        positiveCI = (phat + z*z/(2*total) + z * sqrt((phat*(1-phat)+z*z/(4*total))/total))/(1+z*z/total)
        negativeCI =  (phat + z*z/(2*total) - z * sqrt((phat*(1-phat)+z*z/(4*total))/total))/(1+z*z/total)
        
        return  (phat, positiveCI , negativeCI )


def CountMutations(o, f):
    depths = []

    Aseq = 0
    AtoT = 0
    AtoC = 0
    AtoG = 0

    Tseq = 0
    TtoA = 0
    TtoC = 0
    TtoG = 0

    Cseq = 0
    CtoA = 0
    CtoT = 0
    CtoG = 0

    Gseq = 0
    GtoA = 0
    GtoT = 0
    GtoC = 0

    ins = {}

    dels = {}


    for line in f:
          linebins = line.split()

    #convert sequence information to uppercase
          linebins[4] = linebins[4].replace('t','T')
          linebins[4] = linebins[4].replace('c','C')
          linebins[4] = linebins[4].replace('g','G')
          linebins[4] = linebins[4].replace('a','A')
          linebins[4] = linebins[4].replace('n','N')

    #count depth
          depth = int(linebins[3]) - linebins[4].count('N')

    #skip sites that fall outside of specified start/end ranges, that have insufficient depth or that have clonal mutations or excess frequency of N's:
          if o.end !=0 and int(linebins[1]) < o.start:
                pass
          elif o.end !=0 and int(linebins[1]) > o.end:
                pass
          elif (float(linebins[4].count('N'))/(float(depth) + float(linebins[4].count('N')))) > o.n_cutoff:
                pass
          elif depth < o.mindepth:
                pass
          elif (float(max(linebins[4].count('T'),linebins[4].count('C'),linebins[4].count('G'),linebins[4].count('A'),linebins[4].count('+1'),linebins[4].count('+2'),linebins[4].count('+3'),linebins[4].count('+4'),linebins[4].count('-1'),linebins[4].count('-2'),linebins[4].count('-3'),linebins[4].count('-4'))) / float(depth)) > o.max_clonality:
                pass
          elif (float(max(linebins[4].count('T'),linebins[4].count('C'),linebins[4].count('G'),linebins[4].count('A'),linebins[4].count('+1'),linebins[4].count('+2'),linebins[4].count('+3'),linebins[4].count('+4'),linebins[4].count('-1'),linebins[4].count('-2'),linebins[4].count('-3'),linebins[4].count('-4'))) / float(depth)) < o.min_clonality:
                pass
          else:
            
              #remove N entries
                linebins[4] = linebins[4].replace('N','')
              #remove start line and end line markers
                linebins[4] = re.sub('\$','',linebins[4])
                linebins[4] = re.sub('\^.','',linebins[4])
        
    #count and remove insertions
                newIns = map(int, re.findall(r'\+\d+', linebins[4]))
                if o.unique:
                    newIns = list(set(newIns))
                for length in newIns:
                    if length not in ins:
                        ins[length] = 1
                    else:
                        ins[length] += 1
                    rmStr = r'\+' + str(length) + "."*length
                    linebins[4] = re.sub(rmStr, '', linebins[4])
          
    #count and remove deletions
                newDels = map(str, re.findall(r'-\d+', linebins[4]))
                if o.unique:
                    newDels = list(set(newDels))
                for length in newDels:
                    length = int(length[1:])
                    if length not in dels:
                        dels[length] = 1
                    else:
                        dels[length] += 1
                    rmStr = r'-' + str(length) + "."*length
                    linebins[4] = re.sub(rmStr, '', linebins[4])
         
    #count point mutations
           
                if linebins[2] == 'A':
                      Aseq += depth
                      if linebins[4].count('T') > 0: AtoT += (1 if o.unique else linebins[4].count('T'))
                      if linebins[4].count('C') > 0: AtoC += (1 if o.unique else linebins[4].count('C'))
                      if linebins[4].count('G') > 0: AtoG += (1 if o.unique else linebins[4].count('G'))

                elif linebins[2] == 'T':
                      Tseq += depth
                      if linebins[4].count('A') > 0: TtoA += (1 if o.unique else linebins[4].count('A'))
                      if linebins[4].count('C') > 0: TtoC += (1 if o.unique else linebins[4].count('C'))
                      if linebins[4].count('G') > 0: TtoG += (1 if o.unique else linebins[4].count('G'))
                
                elif linebins[2] == 'C':
                      Cseq += depth
                      if linebins[4].count('A') > 0: CtoA += (1 if o.unique else linebins[4].count('A'))
                      if linebins[4].count('T') > 0: CtoT += (1 if o.unique else linebins[4].count('T'))
                      if linebins[4].count('G') > 0: CtoG += (1 if o.unique else linebins[4].count('G'))
                
                elif linebins[2] == 'G':
                      Gseq += depth
                      if linebins[4].count('A') > 0: GtoA += (1 if o.unique else linebins[4].count('A'))
                      if linebins[4].count('T') > 0: GtoT += (1 if o.unique else linebins[4].count('T'))
                      if linebins[4].count('C') > 0: GtoC += (1 if o.unique else linebins[4].count('C'))

    totalseq = Aseq + Tseq + Cseq + Gseq
    totalptmut = AtoT + AtoC + AtoG + TtoA + TtoC + TtoG + CtoA + CtoT + CtoG + GtoA + GtoT + GtoC
    totalindel = sum(ins) + sum(dels)
                                                
    totalins = sum(ins[n] for n in ins.keys())
    totaldels = sum(dels[n] for n in dels.keys())
                                                
    print ""
    print "Minimum depth", o.mindepth
    print "Maximum clonality", o.max_clonality
    print "Minimum clonality", o.min_clonality
    if o.end != 0:
        print('\nPosition: %s - %s' % (o.start, o.end))
    if o.unique:
        print('\nUnique Counts')
    print ""
    print "A's sequenced", Aseq
    print "A to T:\t", AtoT, '\t%.2e\t%.2e\t%.2e' % Wilson(AtoT,  max(Aseq, 1)) #Output is in the form: Mutation type, number of times mutation is oberseved, frequency, 95% positive CI, 95% negative CI (Confidence Intervals are based on the Wilson Confidence Interval)
    print "A to C:\t", AtoC, '\t%.2e\t%.2e\t%.2e' % Wilson(AtoC,  max(Aseq, 1))
    print "A to G:\t", AtoG, '\t%.2e\t%.2e\t%.2e' % Wilson(AtoG,  max(Aseq, 1))
    print ""
    print "T's sequenced", Tseq
    print "T to A:\t", TtoA, '\t%.2e\t%.2e\t%.2e' % Wilson(TtoA,   max(Tseq, 1))
    print "T to C:\t", TtoC, '\t%.2e\t%.2e\t%.2e' % Wilson(TtoC,  max(Tseq, 1))
    print "T to G:\t", TtoG, '\t%.2e\t%.2e\t%.2e' % Wilson(TtoG,   max(Tseq, 1))
    print ""
    print "C's sequenced", Cseq
    print "C to A:\t", CtoA, '\t%.2e\t%.2e\t%.2e' % Wilson(CtoA,  max(Cseq, 1))
    print "C to T:\t", CtoT, '\t%.2e\t%.2e\t%.2e' % Wilson(CtoT,   max(Cseq, 1))
    print "C to G:\t", CtoG, '\t%.2e\t%.2e\t%.2e' % Wilson(CtoG,   max(Cseq, 1))
    print ""
    print "G's sequenced", Gseq
    print "G to A:\t", GtoA, '\t%.2e\t%.2e\t%.2e' % Wilson(GtoA,   max(Gseq, 1))
    print "G to T:\t", GtoT, '\t%.2e\t%.2e\t%.2e' % Wilson(GtoT,   max(Gseq, 1))
    print "G to C:\t", GtoC, '\t%.2e\t%.2e\t%.2e' % Wilson(GtoC,   max(Gseq, 1))
    print ""
    print "Total nucleotides sequenced", totalseq
    print "Total point mutations", totalptmut
    print "Overall point mutation frequency", '\t%.2e\t%.2e\t%.2e' % Wilson(totalptmut, max(totalseq, 1))
    print ""
    
    insKeys = sorted(ins.items(), key=lambda x: x[0])
    for n in insKeys:
        print('+' + str(n[0]) + ' insertions: ' + str(n[1]) + '\t%.2e\t%.2e\t%.2e' % Wilson(n[1], max(totalseq,1)))
    print ""
    delsKeys = sorted(dels.items(), key=lambda x: x[0])
    for n in delsKeys:
        print('-' + str(n[0]) + ' deletions: ' + str(n[1]) + '\t%.2e\t%.2e\t%.2e' % Wilson(n[1], max(totalseq,1)))
    print ""    
    print "Total insertion events:", totalins
    print "Overall insert frequency:\t%.2e\t%.2e\t%.2e" % Wilson(totalins, max(totalseq, 1))
    print ""
    print "Total deletion events:", totaldels
    print "Overall insert frequency:\t%.2e\t%.2e\t%.2e" % Wilson(totaldels, max(totalseq, 1))

def main():
    parser = ArgumentParser()

    parser.add_argument("-d", "--depth", action="store", type=int, dest="mindepth", 
                      help="Minimum depth for counting mutations at a site (default = 20)", default=20)
    parser.add_argument("-C", "--max_clonality", action="store", type=float, dest="max_clonality",
                      help="Cutoff of mutant reads for scoring a clonal mutation (default = 0.3)", default=0.3)
    parser.add_argument("-c", "--min_clonality", action="store", type=float, dest="min_clonality",
                      help="Cutoff of mutant reads for scoring a clonal mutation (default = 0)", default=0)
    parser.add_argument("-n", "--n_cutoff", action="store", type=float, dest="n_cutoff",
                      help="Maximum fraction of N's allowed to score a position (default = 0.5)", default=0.05)
    parser.add_argument("-s", "--start", action="store", type=int, dest="start",
                      help="Position at which to start scoring for mutations (default = 0)", default=0)
    parser.add_argument("-e", "--end", action="store", type=int, dest="end",
                      help="Position at which to stop scoring for mutations. If set to 0, no position filtering will be performed (default = 0)", default=0)
    parser.add_argument('-u', '--unique', action='store_true', dest='unique', help='run countMutsUnique instead of countMuts')

    o = parser.parse_args()
    f = sys.stdin
    
    CountMutations(o, f)


if __name__ == "__main__":
    main()
