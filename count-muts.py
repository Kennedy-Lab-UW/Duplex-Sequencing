#!/bin/env python

#by Mike Schmitt
#Version 1.3
#January 24, 2013

"""
    
This script pulls out the mutation frequencies from a pileup file given as stdin.

Sites with < 20x depth, and sites with clonal mutations (defined as >30% reads mutated), are excluded from analysis by default.

Usage:

cat seq.pileup | python count-muts.py

"""

from optparse import OptionParser
from math import sqrt
import sys
import re

parser = OptionParser()
parser.add_option("-d", "--depth", action="store", type='int', dest="mindepth", 
                  help="Minimum depth for counting mutations at a site (default = %default)", default=20)
parser.add_option("-c", "--max_clonality", action="store", type='float', dest="max_clonality",
                  help="Upper cutoff of mutant reads for scoring a clonal mutation (default = %default)", default=0.3)
parser.add_option("-C", "--min_clonality", action="store", type='float', dest="min_clonality",
                  help="Lower cutoff of mutant reads for scoring a clonal mutation (default = %default)", default=0)

(options, args) = parser.parse_args()

def Wilson(positive,  total) :
    
        if total == 0:
            
            return 0
            
        freq = float(positive)/float(total)
        z = 1.96 #1.96 = 95%
        phat = float(positive) / total
        positiveCI = (phat + z*z/(2*total) + z * sqrt((phat*(1-phat)+z*z/(4*total))/total))/(1+z*z/total)
        negativeCI =  (phat + z*z/(2*total) - z * sqrt((phat*(1-phat)+z*z/(4*total))/total))/(1+z*z/total)
        
        return  (phat, positiveCI , negativeCI )

f = sys.stdin

#lines = f.readlines()

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

ins1 = 0
ins2 = 0
ins3 = 0
ins4 = 0

del1 = 0
del2 = 0
del3 = 0
del4 = 0


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
    
#skip sites that have insufficient depth or that have clonal mutations    
      if depth < options.mindepth:
            pass
      elif (float(max(linebins[4].count('T'),linebins[4].count('C'),linebins[4].count('G'),linebins[4].count('A'),linebins[4].count('+1'),linebins[4].count('+2'),linebins[4].count('+3'),linebins[4].count('+4'),linebins[4].count('-1'),linebins[4].count('-2'),linebins[4].count('-3'),linebins[4].count('-4'))) / float(depth)) > options.max_clonality:
            pass
      elif (float(max(linebins[4].count('T'),linebins[4].count('C'),linebins[4].count('G'),linebins[4].count('A'),linebins[4].count('+1'),linebins[4].count('+2'),linebins[4].count('+3'),linebins[4].count('+4'),linebins[4].count('-1'),linebins[4].count('-2'),linebins[4].count('-3'),linebins[4].count('-4'))) / float(depth)) < options.min_clonality:
            pass
      else:          
        
          #remove start line, end line, and N entries
            linebins[4] = re.sub('\$','',linebins[4])
            linebins[4] = re.sub('\^.','',linebins[4])    
            linebins[4] = linebins[4].replace('N','')          
          
#count and remove insertions
            ins1 += linebins[4].count('+1')
            ins2 += linebins[4].count('+2')
            ins3 += linebins[4].count('+3')
            ins4 += linebins[4].count('+4')

            linebins[4] = re.sub('\+1.','',linebins[4])    
            linebins[4] = re.sub('\+2..','',linebins[4])
            linebins[4] = re.sub('\+3...','',linebins[4])    
            linebins[4] = re.sub('\+4....','',linebins[4])    
      
#count and remove deletions
            del1 += linebins[4].count('-1')
            del2 += linebins[4].count('-2')
            del3 += linebins[4].count('-3')
            del4 += linebins[4].count('-4')

            linebins[4] = linebins[4].replace('*','')    
            linebins[4] = re.sub('-1.','',linebins[4])    
            linebins[4] = re.sub('-2..','',linebins[4])    
            linebins[4] = re.sub('-3...','',linebins[4])    
            linebins[4] = re.sub('-4....','',linebins[4])    
     
#count point mutations
       
            if linebins[2] == 'A':
                  Aseq += depth
                  AtoT += linebins[4].count('T')
                  AtoC += linebins[4].count('C')
                  AtoG += linebins[4].count('G')

            elif linebins[2] == 'T':
                  Tseq += depth
                  TtoA += linebins[4].count('A')
                  TtoC += linebins[4].count('C')
                  TtoG += linebins[4].count('G')
            
            elif linebins[2] == 'C':
                  Cseq += depth
                  CtoA += linebins[4].count('A')
                  CtoT += linebins[4].count('T')
                  CtoG += linebins[4].count('G')
            
            elif linebins[2] == 'G':
                  Gseq += depth
                  GtoA += linebins[4].count('A')
                  GtoT += linebins[4].count('T')
                  GtoC += linebins[4].count('C')

totalseq = Aseq + Tseq + Cseq + Gseq
totalptmut = AtoT + AtoC + AtoG + TtoA + TtoC + TtoG + CtoA + CtoT + CtoG + GtoA + GtoT + GtoC
totalins = ins1 + ins2 + ins3 + ins4
totaldels = del1 + del2 + del3 + del4                                            
  
print ""
print "Minimum depth", options.mindepth
print "Maximum clonality", options.max_clonality
print "Minimum clonality", options.min_clonality
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
print ""
print "Total point mutations", totalptmut
print ""
print "Overall point mutation frequency", '\t%.2e\t%.2e\t%.2e' % Wilson(totalptmut, max(totalseq, 1))
print ""
print ""
print "+1 insertions", ins1, '\t%.2e\t%.2e\t%.2e' % Wilson(ins1, max(totalseq,1))
print "+2 insertions", ins2, '\t%.2e\t%.2e\t%.2e' % Wilson(ins2, max(totalseq,1))
print "+3 insertions", ins3, '\t%.2e\t%.2e\t%.2e' % Wilson(ins3, max(totalseq,1))
print "+4 insertions", ins4, '\t%.2e\t%.2e\t%.2e' % Wilson(ins4, max(totalseq,1))
print ""
print "-1 deletions",del1, '\t%.2e\t%.2e\t%.2e' % Wilson(del1, max(totalseq,1))
print "-2 deletions",del2, '\t%.2e\t%.2e\t%.2e' % Wilson(del2, max(totalseq,1))
print "-3 deletions",del3, '\t%.2e\t%.2e\t%.2e' % Wilson(del3, max(totalseq,1))
print "-4 deletions ",del4, '\t%.2e\t%.2e\t%.2e' % Wilson(del4, max(totalseq,1))
print ""
print "Total insertion events:", totalins
print "Overall insert frequency:\t%.2e\t%.2e\t%.2e" % Wilson(totalins, max(totalseq, 1))
print ""
print "Total deletion events:", totaldels
print "Overall insert frequency:\t%.2e\t%.2e\t%.2e" % Wilson(totaldels, max(totalseq, 1))




'''
totalseq = Aseq + Tseq + Cseq + Gseq
totalptmut = AtoT + AtoC + AtoG + TtoA + TtoC + TtoG + CtoA + CtoT + CtoG + GtoA + GtoT + GtoC
totalindel = ins1 + ins2 + ins3 + ins4 + del1 + del2 + del3 + del4
                                            
ptmutfreq = (float (totalptmut) / float(totalseq))
indelfreq = (float (totalindel) / float(totalseq))
                                            
print ""
print "Minimum depth", options.mindepth
print "Maximum clonality", options.max_clonality
print "Minimum clonality", options.min_clonality
print ""
print "A's sequenced", Aseq
print "A to T:", AtoT, ' %.2e' % (float(AtoT) / float(max(Aseq,1)))
print "A to C:", AtoC, ' %.2e' % (float(AtoC) / float(max(Aseq,1)))
print "A to G:", AtoG, ' %.2e' % (float(AtoG) / float(max(Aseq,1)))
print ""
print "T's sequenced", Tseq
print "T to A:", TtoA, ' %.2e' % (float(TtoA) / float(max(Tseq,1)))
print "T to C:", TtoC, ' %.2e' % (float(TtoC) / float(max(Tseq,1)))
print "T to G:", TtoG, ' %.2e' % (float(TtoG) / float(max(Tseq,1)))
print ""
print "C's sequenced", Cseq
print "C to A:", CtoA, ' %.2e' % (float(CtoA) / float(max(Cseq,1)))
print "C to T:", CtoT, ' %.2e' % (float(CtoT) / float(max(Cseq,1)))
print "C to G:", CtoG, ' %.2e' % (float(CtoG) / float(max(Cseq,1)))
print ""
print "G's sequenced", Gseq
print "G to A:", GtoA, ' %.2e' % (float(GtoA) / float(max(Gseq,1)))
print "G to T:", GtoT, ' %.2e' % (float(GtoT) / float(max(Gseq,1)))
print "G to C:", GtoC, ' %.2e' % (float(GtoC) / float(max(Gseq,1)))
print ""
print "+1 insertions", ins1, ' %.2e' % (float(ins1) / float(max(totalseq,1)))
print "-1 deletions",del1, ' %.2e' % (float(del1) / float(max(totalseq,1)))
print "Total indel events", totalindel
print "Overall indel frequency", '%.2e' % indelfreq
print ""
print "Total nucleotides sequenced", totalseq
print "Total point mutations", totalptmut
print "Overall point mutation frequency", '%.2e' % ptmutfreq
'''
