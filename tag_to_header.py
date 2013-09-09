#By Joe Hiatt, Scott Kennedy and Mike Schmitt
#Version 1.21
#February 5, 2013

import sys
from optparse import OptionParser

opts = OptionParser()
opts.add_option('', '--adapter1',  type='string',  default=None,  dest='adapterSeq1')
opts.add_option('', '--adapter2',  type='string',  default=None,  dest='adapterSeq2')
opts.add_option('','--infile1',type='string',default=None,dest='infile1')
opts.add_option('','--infile2',type='string',default=None,dest='infile2')
opts.add_option('','--outfile1',type='string',default=None,dest='outfile1')
opts.add_option('','--outfile2',type='string',default=None,dest='outfile2')
o,args=opts.parse_args()

inFile1 = open(o.infile1,  'r')
inFile2 = open(o.infile2,  'r')
outFile1 = open(o.outfile1,  'w')
outFile2 = open(o.outfile2,  'w')

 #this is the function that extracts the UID tags from both the forward and reverse read.  Assigns read1 the sequence from some position to the end, then read2 from some position to the end
 #then assigns tag1 from the 5'-end to length of the UID tag for read1 and then read 2.  Change as necessary to your appropriate adapter design.
tagExtractFxn ='lambda x:(x[0][17:], x[1][17:], x[0][:12], x[1][:12])'

#this function renames the header with the formatting of *header coordinates,etc*, *tag from read1*, *tag from read2*, *read designation from original header*
hdrRenameFxn ='lambda x,y,z:"%s#%s%s/%s" % (x.split("#")[0], y, z, x.split("/")[-1])'

tagExtractFxn = eval (tagExtractFxn)
hdrRenameFxn = eval (hdrRenameFxn)

infilelist = [inFile1,  inFile2] 
outfilelist = [ outFile1, outFile2 ]

ctr = 0
totalreads = 0
goodreads = 0
badtag = 0
nospacer = 0

for line1 in infilelist[0]:
    line1 = line1.strip()
    line2 = infilelist[1].next().strip()

    ctr += 1
    
    if ctr%4==1:
        hdr1 = line1
        hdr2 = line2

    elif ctr%4==2:
        read1 = line1
        read2 = line2

    elif ctr%4==0:
        totalreads += 1

        if o.adapterSeq1 != None and o.adapterSeq2 != None and ( read1[12:16] != o.adapterSeq1 or read2[12:16] != o.adapterSeq2 )  : #filtering logic is contained here changes the values as necessary to the location of your adapter spacers (0 indexed)
            nospacer += 1

        else :
            #extract tags
            read1, read2,  tag1,  tag2 = tagExtractFxn((read1, read2))
            line1, line2,  trash,  trash2 = tagExtractFxn (( line1, line2 ))
            
            #header reconstruction
            hdr1 = hdrRenameFxn ( hdr1, tag1, tag2 ) 
            hdr2 = hdrRenameFxn ( hdr2, tag1, tag2 )
            
            #tag1Hash = int(tag1.count('A')) + int(tag1.count('C')) + int(tag1.count('G')) + int(tag1.count('T'))
            #tag2Hash = int(tag1.count('A')) + int(tag1.count('C')) + int(tag1.count('G')) + int(tag1.count('T'))
            
            #fastq reconstruction
            if (tag1.isalpha() and tag1.count('N') == 0) and (tag2.isalpha() and tag2.count('N') == 0) :
                outfilelist[ 0 ].write( "%s\n%s\n+\n%s\n"%( hdr1, read1[4:], line1[4:] ) ) # keeps nucleotides 5 onwards (1-indexed) from read that is already trimmed of UID and adapter spacer
                outfilelist[ 1 ].write( "%s\n%s\n+\n%s\n"%( hdr2, read2[4:], line2[4:] ) ) # keeps nucleotides 5 onwards (1-indexed) from read that is already trimmed of UID and adapter spacer
                #outfilelist[ 0 ].write( "%s\n%s\n+\n%s\n"%( hdr1, read1, line1 ) ) # keeps nucleotides 8 onwards (1-indexed) from read that is already trimmed of UID and adapter spacer
                #outfilelist[ 1 ].write( "%s\n%s\n+\n%s\n"%( hdr2, read2, line2 ) ) # keeps nucleotides 8 onwards (1-indexed) from read that is already trimmed of UID and adapter spacer

                goodreads += 1

            else: badtag += 1
    
    if ctr%400000==0:
        sys.stderr.write( "processed %d reads...\n"%(ctr/4) )
        print "Total sequences processed:", totalreads
        print "Sequences passing filter:", goodreads
        print "Missing spacers:", nospacer
        print "Bad tags:", badtag
        continue
    
for file in infilelist:
    file.close()

for file in outfilelist:
    file.close()

print " "
print "Summary statistics:"
print "Total sequences processed:", totalreads
print "Sequences passing filter:", goodreads
print "Missing spacers:", nospacer
print "Bad tags:", badtag
