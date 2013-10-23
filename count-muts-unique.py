#!/bin/env python

#by Mike Schmitt
#Version 1.3
#Jan 31, 2013
#Modified from count-muts.py v1.1
#Edited by Brendan Kohrn to fix a problem with 0 depth where no 0 depth should be

"""
    
This script pulls out the mutation frequencies from a pileup file given as stdin.

Sites with < 20x depth, and sites with clonal mutations (defined as >30% reads mutated), are excluded from analysis by default.

This version counts each mutation exactly once (i.e. clonal expansions are counted as a single mutation)

Usage:

cat seq.pileup | python count-muts.py

"""

from optparse import OptionParser
import sys
import re

parser = OptionParser()

parser.add_option("-d", "--depth", action="store", type='int', dest="mindepth",
                  help="Minimum depth for counting mutations at a site (default = %default)", default=20)
parser.add_option("-c", "--max_clonality", action="store", type='float', dest="max_clonality",
                  help="Cutoff of mutant reads for scoring a clonal mutation (default = %default)", default=0.3)
parser.add_option("-C", "--min_clonality", action="store", type='float', dest="min_clonality",
                  help="Cutoff of mutant reads for scoring a clonal mutation (default = %default)", default=0)
parser.add_option("-n", "--n_cutoff", action="store", type='float', dest="n_cutoff",
                  help="Maximum fraction of N's allowed to score a position (default = %default)", default=0.05)
parser.add_option("-s", "--start", action="store", type='int', dest="start",
                  help="Position at which to start scoring for mutations (default = $default)", default=0)
parser.add_option("-e", "--end", action="store", type='int', dest="end",
                  help="Position at which to stop scoring for mutations. If set to 0, no position filtering will be performed (default = $default)", default=0)

(options, args) = parser.parse_args()

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

#remove start line and end line markers
      linebins[4] = re.sub('\$','',linebins[4])
      linebins[4] = re.sub('\^.','',linebins[4])

#count depth
      depth = int(linebins[3]) - linebins[4].count('N')

#skip sites that fall outside of specified start/end ranges, that have insufficient depth or that have clonal mutations or excess frequency of N's:
      if options.end !=0 and int(linebins[1]) < options.start:
            pass
      elif options.end !=0 and int(linebins[1]) > options.end:
            pass
      elif (float(linebins[4].count('N'))/(float(depth) + float(linebins[4].count('N')))) > options.n_cutoff:
            pass
      elif depth < options.mindepth:
            pass
      elif (float(max(linebins[4].count('T'),linebins[4].count('C'),linebins[4].count('G'),linebins[4].count('A'),linebins[4].count('+1'),linebins[4].count('+2'),linebins[4].count('+3'),linebins[4].count('+4'),linebins[4].count('-1'),linebins[4].count('-2'),linebins[4].count('-3'),linebins[4].count('-4'))) / float(depth)) > options.max_clonality:
            pass
      elif (float(max(linebins[4].count('T'),linebins[4].count('C'),linebins[4].count('G'),linebins[4].count('A'),linebins[4].count('+1'),linebins[4].count('+2'),linebins[4].count('+3'),linebins[4].count('+4'),linebins[4].count('-1'),linebins[4].count('-2'),linebins[4].count('-3'),linebins[4].count('-4'))) / float(depth)) < options.min_clonality:
            pass
      else:          
        
          #remove N entries
            linebins[4] = linebins[4].replace('N','')          
          
#count and remove insertions
            if linebins[4].count('+1') > 0: ins1 += 1
            if linebins[4].count('+2') > 0: ins1 += 1
            if linebins[4].count('+3') > 0: ins1 += 1
            if linebins[4].count('+4') > 0: ins1 += 1

            linebins[4] = re.sub('\+1.','',linebins[4])    
            linebins[4] = re.sub('\+2..','',linebins[4])
            linebins[4] = re.sub('\+3...','',linebins[4])    
            linebins[4] = re.sub('\+4....','',linebins[4])    
      
#count and remove deletions
            if linebins[4].count('-1'): del1 += 1
            if linebins[4].count('-2'): del2 += 1
            if linebins[4].count('-3'): del3 += 1
            if linebins[4].count('-4'): del4 += 1

            linebins[4] = linebins[4].replace('*','')    
            linebins[4] = re.sub('-1.','',linebins[4])    
            linebins[4] = re.sub('-2..','',linebins[4])    
            linebins[4] = re.sub('-3...','',linebins[4])    
            linebins[4] = re.sub('-4....','',linebins[4])    
     
#count point mutations
       
            if linebins[2] == 'A':
                  Aseq += depth
                  if linebins[4].count('T') > 0: AtoT += 1
                  if linebins[4].count('C') > 0: AtoC += 1
                  if linebins[4].count('G') > 0: AtoG += 1

            elif linebins[2] == 'T':
                  Tseq += depth
                  if linebins[4].count('A') > 0: TtoA += 1
                  if linebins[4].count('C') > 0: TtoC += 1
                  if linebins[4].count('G') > 0: TtoG += 1
            
            elif linebins[2] == 'C':
                  Cseq += depth
                  if linebins[4].count('A') > 0: CtoA += 1
                  if linebins[4].count('T') > 0: CtoT += 1
                  if linebins[4].count('G') > 0: CtoG += 1
            
            elif linebins[2] == 'G':
                  Gseq += depth
                  if linebins[4].count('A') > 0: GtoA += 1
                  if linebins[4].count('T') > 0: GtoT += 1
                  if linebins[4].count('C') > 0: GtoC += 1

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
