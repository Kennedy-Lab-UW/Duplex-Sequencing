"""
Count Muts Unique

by Brendan Kohrn and Mike Schmitt
Version 1.41
August 20, 2014
Modified from count-muts.py and count-muts-unique.py
Edited by Brendan Kohrn to fix a problem with 0 depth where no 0 depth should be and to allow n-length indels, and to merge the two conditions into one program.  

This script pulls out the mutation frequencies from a pileup file given as stdin, or can take an imput file using the -i option, and writes to stdout, or can take an output file name.   

Sites with less than MINDEPTH, or clonalities outside of the range MIN_CLONALITY-MAX_CLONALITY, are excluded from analysis.

If -u is specified, this program counts each mutation exactly once (i.e. clonal expansions are counted as a single mutation)

Usage:

cat seq.pileup | CountMuts.py [-h] [-d MINDEPTH] [-C MAX_CLONALITY] [-c MIN_CLONALITY] [-n N_CUTOFF] [-s START] [-e END] [-u] > outfile.countmuts

optional arguments:
  -h, --help            show this help message and exit
  -d MINDEPTH, --depth MINDEPTH
                        Minimum depth for counting mutations at a site
                        (default = 20)
  -C MAX_CLONALITY, --max_clonality MAX_CLONALITY
                        Cutoff of mutant reads for scoring a clonal mutation
                        (default = 0.3)
  -c MIN_CLONALITY, --min_clonality MIN_CLONALITY
                        Cutoff of mutant reads for scoring a clonal mutation
                        (default = 0)
  -n N_CUTOFF, --n_cutoff N_CUTOFF
                        Maximum fraction of N's allowed to score a position
                        (default = 0.05)
  -s START, --start START
                        Position at which to start scoring for mutations
                        (default = 0)
  -e END, --end END     Position at which to stop scoring for mutations. If
                        set to 0, no position filtering will be performed
                        (default = 0)
  -u, --unique          run countMutsUnique instead of countMuts

"""

from __future__ import print_function
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


def CountMutations(o, f, fOut):
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

    ins = {0:0}

    dels = {0:0}

    #mpFile = open("testMPfile.mutpos", "w") #ADDED
    #mpFirst = True #ADDED
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

    #skip sites that fall outside of specified start/end ranges, that have insufficient depth or that have clonal mutations or excess frequency of N's:
          if o.end !=0 and int(linebins[1]) < o.start:
                pass
          elif o.end !=0 and int(linebins[1]) > o.end:
                pass
          elif (float(linebins[4].count('N'))/(float(depth) + float(linebins[4].count('N')))) > o.n_cutoff:
                pass
          elif depth < o.mindepth:
                pass
          elif (float(max(linebins[4].count('T'),linebins[4].count('C'),linebins[4].count('G'),linebins[4].count('A'), (max(newIns.count(n) for n in list(set(newIns))) if newIns != [] else 0), (max(newDels.count(m) for m in list(set(newDels))) if newDels != [] else 0))) / float(depth)) > o.max_clonality:
                pass
          elif (float(max(linebins[4].count('T'),linebins[4].count('C'),linebins[4].count('G'),linebins[4].count('A'), (max(newIns.count(n) for n in list(set(newIns))) if newIns != [] else 0), (max(newDels.count(m) for m in list(set(newDels))) if newDels != [] else 0))) / float(depth)) < o.min_clonality:
                pass
          else:
            
              #remove N entries
                #linebins[4] = linebins[4].replace('N','')
              #remove start line and end line markers
                linebins[4] = re.sub('\$','',linebins[4])
                linebins[4] = re.sub('\^.','',linebins[4])
         
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
                #if mpFirst: #ADDED
                    #mpFirst = False #ADDED
                #else: #ADDED
                    #mpFile.write("\n") #ADDED
                #mpFile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (linebins[0],linebins[2], linebins[1], depth, linebins[4].count('T') + linebins[4].count('C') + linebins[4].count('G') + linebins[4].count('A'), linebins[4].count('T'), linebins[4].count('C'), linebins[4].count('G'), linebins[4].count('A'), len(newIns), len(newDels), linebins[4].count('N'))) #ADDED

    totalseq = Aseq + Tseq + Cseq + Gseq
    totalptmut = AtoT + AtoC + AtoG + TtoA + TtoC + TtoG + CtoA + CtoT + CtoG + GtoA + GtoT + GtoC
    totalindel = sum(ins) + sum(dels)
                                                
    totalins = sum(ins[n] for n in ins.keys())
    totaldels = sum(dels[n] for n in dels.keys())
    #mpFile.close() #ADDED

    print("\nMinimum depth: %s" % o.mindepth, file = fOut)
    print("Clonality: %s - %s" % (o.min_clonality, o.max_clonality), file = fOut)
    if o.end != 0: 
        print('Position: %s - %s' % (o.start, o.end), file = fOut)
    if o.unique:
        print('Unique Counts', file = fOut)
    print("\nA's sequenced: %s" % Aseq, file = fOut)
    print(("A to T:\t%s" % AtoT) + ('\t%.2e\t%.2e\t%.2e' % Wilson(AtoT,  max(Aseq, 1))), file = fOut) #Output is in the form: Mutation type, number of times mutation is oberseved, frequency, 95% positive CI, 95% negative CI (Confidence Intervals are based on the Wilson Confidence Interval)
    print(("A to C:\t%s" % AtoC) + ('\t%.2e\t%.2e\t%.2e' % Wilson(AtoC,  max(Aseq, 1))), file = fOut)
    print(("A to G:\t%s" % AtoG) + ('\t%.2e\t%.2e\t%.2e' % Wilson(AtoG,  max(Aseq, 1))), file = fOut)
    print("\nT's sequenced: %s" % Tseq, file = fOut)
    print(("T to A:\t%s" % TtoA) + ('\t%.2e\t%.2e\t%.2e' % Wilson(TtoA,   max(Tseq, 1))), file = fOut)
    print(("T to C:\t%s" % TtoC) + ('\t%.2e\t%.2e\t%.2e' % Wilson(TtoC,  max(Tseq, 1))), file = fOut)
    print(("T to G:\t%s" % TtoG) + ('\t%.2e\t%.2e\t%.2e' % Wilson(TtoG,   max(Tseq, 1))), file = fOut)
    print("\nC's sequenced: %s" % Cseq, file = fOut)
    print(("C to A:\t%s" % CtoA) + ('\t%.2e\t%.2e\t%.2e' % Wilson(CtoA,  max(Cseq, 1))), file = fOut)
    print(("C to T:\t%s" % CtoT) + ('\t%.2e\t%.2e\t%.2e' % Wilson(CtoT,   max(Cseq, 1))), file = fOut)
    print(("C to G:\t%s" % CtoG) + ('\t%.2e\t%.2e\t%.2e' % Wilson(CtoG,   max(Cseq, 1))), file = fOut)
    print("\nG's sequenced: %s" % Gseq, file = fOut)
    print(("G to A:\t%s" % GtoA) + ('\t%.2e\t%.2e\t%.2e' % Wilson(GtoA,   max(Gseq, 1))), file = fOut)
    print(("G to T:\t%s" % GtoT) + ('\t%.2e\t%.2e\t%.2e' % Wilson(GtoT,   max(Gseq, 1))), file = fOut)
    print(("G to C:\t%s" % GtoC) + ('\t%.2e\t%.2e\t%.2e' % Wilson(GtoC,   max(Gseq, 1))), file = fOut)
    print("\nTotal nucleotides sequenced: %s" % totalseq, file = fOut)
    print("Total point mutations: %s" % totalptmut, file = fOut)
    print('Overall point mutation frequency:\t%.2e\t%.2e\t%.2e\n' % Wilson(totalptmut, max(totalseq, 1)), file = fOut)
    
    insKeys = sorted(ins.items(), key=lambda x: x[0])
    for n in insKeys:
        print(('+%s insertions: %s' % (n[0], n[1])) + ('\t%.2e\t%.2e\t%.2e' % Wilson(n[1], max(totalseq,1))), file = fOut)
    if dels != {}:
        print('', file = fOut)
    delsKeys = sorted(dels.items(), key=lambda x: x[0])
    for n in delsKeys:
        print(('-%s deletions: %s' % (n[0], n[1])) + ('\t%.2e\t%.2e\t%.2e' % Wilson(n[1], max(totalseq,1))), file = fOut)   
    print("\nTotal insertion events: %s" % totalins, file = fOut)
    print("Overall insert frequency:\t%.2e\t%.2e\t%.2e" % Wilson(totalins, max(totalseq, 1)), file = fOut)
    print("\nTotal deletion events: %s" % totaldels, file = fOut)
    print("Overall deletion frequency:\t%.2e\t%.2e\t%.2e" % Wilson(totaldels, max(totalseq, 1)), file = fOut)

def main():
    parser = ArgumentParser()
    parser.add_argument('-i', '--infile', action ='store', dest = 'inFile', help = 'An imput file. If None, defaults to stdin. [None]', default = None)
    parser.add_argument('-o', '--outfile', action = 'store', dest = 'outFile', help = 'A filename for the output file.  If None, outputs to stdout.  [None]', default = None)
    parser.add_argument("-d", "--depth", action="store", type=int, dest="mindepth", 
                      help="Minimum depth for counting mutations at a site [%(default)s]", default=20)
    parser.add_argument("-c", "--min_clonality", action="store", type=float, dest="min_clonality",
                      help="Cutoff of mutant reads for scoring a clonal mutation [%(default)s]", default=0)
    parser.add_argument("-C", "--max_clonality", action="store", type=float, dest="max_clonality",
                      help="Cutoff of mutant reads for scoring a clonal mutation [%(default)s]", default=0.3)
    parser.add_argument("-n", "--n_cutoff", action="store", type=float, dest="n_cutoff",
                      help="Maximum fraction of N's allowed to score a position [%(default)s]", default=0.05)
    parser.add_argument("-s", "--start", action="store", type=int, dest="start",
                      help="Position at which to start scoring for mutations [%(default)s]", default=0)
    parser.add_argument("-e", "--end", action="store", type=int, dest="end",
                      help="Position at which to stop scoring for mutations. If set to 0, no position filtering will be performed [%(default)s]", default=0)
    parser.add_argument('-u', '--unique', action='store_true', dest='unique', help='Run countMutsUnique instead of countMuts')

    o = parser.parse_args()
    if o.inFile != None:
        f = open(o.inFile, 'r')
    else:
        f = sys.stdin
    if o.outFile != None:
        fOut = open(o.outFile, 'w')
    else:
        fOut = sys.stdout
    CountMutations(o, f, fOut)


if __name__ == "__main__":
    main()
