#!/bin/env python

#By Mike Schmitt
#Version 1.1
#June 14, 2012

"""
This script gives position-specific mutation frequencies from a tagcounts file given as stdin.

The output is tab-delimited and specifies:
chromosome number, template base, nucleotide position, depth, mutations to T, C, G, A, insertions, deletions

Usage:

cat seq.pileup | python count-muts.py

Mike Schmitt, Jan 9, 2012
"""

from optparse import OptionParser
import sys
import re
import csv
import string

parser = OptionParser()
parser.add_option("-d", "--depth", action="store", type='int', dest="mindepth", 
                  help="Minimum depth for counting mutations at a site (default = %default)", default=1)
parser.add_option("-n", "--num_muts", action="store", type='int', dest="num_muts", 
                  help="Minimum number of mutations for scoring a site (default = %default)", default=0)
parser.add_option("-C", "--clonality_max", action="store", type='float', dest="clonal_max",
                  help="Maximum clonality to score a mutation (default = %default)", default=1)
parser.add_option("-c", "--clonality_min", action="store", type='float', dest="clonal_min",
                  help="Minimum clonality to score a mutation (default = %default)", default=0)

(options, args) = parser.parse_args()

f = sys.stdin

lines = f.readlines()

chrom=[]
pos=[]
muts = []
depths = []
template =[]
Tcount=[]
Ccount=[]
Gcount=[]
Acount=[]
inscount=[]
delcount=[]
for i,line in enumerate( lines ):

      linebins = line.split()

#convert sequence information to uppercase
      linebins[4] = linebins[4].replace('t','T')
      linebins[4] = linebins[4].replace('c','C')
      linebins[4] = linebins[4].replace('g','G')
      linebins[4] = linebins[4].replace('a','A')
      linebins[4] = linebins[4].replace('n','N')

#remove start line, end line, and N entries, as well as 1st and last nucleotide of a read.
      linebins[4] = re.sub('.\$','',linebins[4])
      linebins[4] = re.sub('\^..','',linebins[4])    
      linebins[4] = linebins[4].replace('N','')      

#count and remove insertions
      ins1 = linebins[4].count('+1')
      ins2 = linebins[4].count('+2')
      ins3 = linebins[4].count('+3')
      ins4 = linebins[4].count('+4')

      linebins[4] = re.sub('\+1.','',linebins[4])    
      linebins[4] = re.sub('\+2..','',linebins[4])
      linebins[4] = re.sub('\+3...','',linebins[4])    
      linebins[4] = re.sub('\+4....','',linebins[4])    
      
#count and remove deletions
      del1 = linebins[4].count('-1')
      del2 = linebins[4].count('-2')
      del3 = linebins[4].count('-3')
      del4 = linebins[4].count('-4')

      linebins[4] = linebins[4].replace('*','')    
      linebins[4] = re.sub('-1.','',linebins[4])    
      linebins[4] = re.sub('-2..','',linebins[4])    
      linebins[4] = re.sub('-3...','',linebins[4])    
      linebins[4] = re.sub('-4....','',linebins[4])    

#count depth                                                                                                
      depth = len(linebins[4])
                                                    
#skip lines that do not meet filtering criteria
      if    (
            depth < options.mindepth
            or
            ((float(max(linebins[4].count('T'),linebins[4].count('C'),linebins[4].count('G'),linebins[4].count('A'),ins1,ins2,ins3,ins4,del1,del2,del3,del4)) / float(depth)) > options.clonal_max)
            or
            ((float(max(linebins[4].count('T'),linebins[4].count('C'),linebins[4].count('G'),linebins[4].count('A'),linebins[4].count('+1'),linebins[4].count('+2'),linebins[4].count('+3'),linebins[4].count('+4'),linebins[4].count('-1'),linebins[4].count('-2'),linebins[4].count('-3'),linebins[4].count('-4'))) / float(depth)) < options.clonal_min)
            or    
            (max(float(linebins[4].count('T')),float(linebins[4].count('C')),float(linebins[4].count('G')),float(linebins[4].count('A')),float(linebins[4].count('+1')),float(linebins[4].count('+2')),float(linebins[4].count('+3')),float(linebins[4].count('+4')),float(linebins[4].count('-1')),float(linebins[4].count('-2')),float(linebins[4].count('-3')),float(linebins[4].count('-4'))) < options.num_muts) 
            ):
            pass
                                                    
      else:

#count position-specific mutation frequency
                                                    
            mut = linebins[4].count('T') + linebins[4].count('C') + linebins[4].count('G') + linebins[4].count('A')
            Tcount.append(linebins[4].count('T'))
            Ccount.append(linebins[4].count('C'))
            Gcount.append(linebins[4].count('G'))
            Acount.append(linebins[4].count('A'))
            inscount.append(ins1+ins2+ins3+ins4)
            delcount.append(del1+del2+del3+del4)
            chrom.append(linebins[0])
            pos.append(linebins[1])
            template.append(linebins[2])
            depths.append(depth)
            muts.append(mut)
                                                  
script_output=zip(chrom, template, pos, depths, muts, Tcount, Ccount, Gcount, Acount, inscount, delcount)

csv_writer = csv.writer(sys.stdout, delimiter='\t')
csv_writer.writerows(script_output)
